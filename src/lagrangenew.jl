"""
  lagrangesolve(graph)
  solve graph using lagrange decomposition
"""
function lagrangesolve(graph;max_iterations=10,update_method=:subgradient,ϵ=0.1,timelimit=3600,α=2,lagrangeheuristic=fixbinaries)
  n = normalizegraph(graph)
  lgprepare(graph)

  starttime = time()
  s = Solution(method=:dual_decomposition)

  nmult = graph.attributes[:numlinks]
  λ = zeros(nmult) # Array{Float64}(nmult)
  x = zeros(nmult,2) # Linking variables values
  res = zeros(nmult) # Residuals
  nodes = [node for node in values(getnodes(graph))]
  graph.attributes[:α] = α

  for iter in 1:max_iterations
    iterstart = time()
    # Solve subproblems
    Zk = 0
    for node in nodes
       Zk += solvenode(node,λ,x)
    end
    graph.attributes[:Zk] = Zk

    # Update residuals
    res = x[:,1] - x[:,2]

    # Check convergence
    if mean(res) < ϵ
      println("Converged")
      break
    end

    # Update multipliers
    iterval = updatemultipliers(graph,λ,res,update_method,lagrangeheuristic)

    # Save summary
    itertime = time() - iterstart
    tstamp = time() - starttime
    saveiteration(s,tstamp,[iterval,Zk,itertime],n)
  end
  return s
end

# Preprocess function
"""
  lgprepare(graph::PlasmoGraph)
  Prepares the graph to apply lagrange decomposition algorithm
"""
function lgprepare(graph::PlasmoGraph)
  if haskey(graph.attributes,:preprocessed)
    return true
  end

  links = getlinkconstraints(graph)
  nmult = length(links) # Number of multipliers
  graph.attributes[:numlinks] = nmult
  graph.attributes[:sense] = getmodel(getnodes(graph)[1]).objSense
  sense = graph.attributes[:sense]
  graph.attributes[:mflat] = create_flat_graph_model(graph)
  graph.attributes[:cuts] = []

  # Create Lagrange Master
  ms = Model(solver=graph.solver)
  @variable(ms, η)
  @variable(ms, λ[1:nmult])
  mssense = :Min
  if sense == :Min
    mssense = :Max
  end
  @objective(ms, mssense, η)

  graph.attributes[:lgmaster] = ms

  # Each node most save its initial objective
  for n in values(getnodes(graph))
    mn = getmodel(n)
    mn.ext[:preobj] = mn.obj
    mn.ext[:multmap] = Dict()
    mn.ext[:varmap] = Dict()
  end

  # Maps
  # Multiplier map to know which component of λ to take
  # Varmap knows what values to post where
  for (i,lc) in enumerate(links)
    for j in 1:length(lc.terms.vars)
      var = lc.terms.vars[j]
      var.m.ext[:multmap][i] = (lc.terms.coeffs[j],lc.terms.vars[j])
      var.m.ext[:varmap][var] = (i,j)
    end
  end

  graph.attributes[:preprocessed] = true
end

# Solve a single subproblem
function solvenode(node,λ,x)
  m = getmodel(node)
  m.obj = m.ext[:preobj]
  # Add dualized part to objective function
  for k in keys(m.ext[:multmap])
    m.obj += λ[k]*m.ext[:multmap][k][1]*m.ext[:multmap][k][2]
  end

  # Optional: If my residuals are zero, do nothing

  solve(m)
  for v in keys(m.ext[:varmap])
    val = getvalue(v)
    x[m.ext[:varmap][v]...] = val
  end

  objval = getobjectivevalue(m)
  node.attributes[:objective] = objval
  node.attributes[:solvetime] = getsolvetime(m)

  return objval
end

function updatemultipliers(graph,λ,res,method,lagrangeheuristic=nothing)
  if method == :subgradient
    subgradient(graph,λ,res,lagrangeheuristic)
  elseif method == :cuttingplanes
    cuttingplanes(graph,λ,res)
  end
end

function subgradient(graph,λ,res,lagrangeheuristic)
  α = graph.attributes[:α]
  bound = lagrangeheuristic(graph)
  Zk = graph.attributes[:Zk]
  λ += α*abs(Zk-bound)/(norm(res)^2)*res
  return bound
end

function cuttingplanes(graph,λ,res)
  ms = graph.attributes[:lgmaster]
  Zk = graph.attributes[:Zk]
  nmult = graph.attributes[:numlinks]

  λvar = getindex(ms, :λ)
  η = getindex(ms,:η)

  sense = ms.objSense
  if sense == :Min
    @constraint(ms, η >= Zk + sum(λvar[j]*res[j] for j in 1:nmult))
  else
    @constraint(ms, η <= Zk + sum(λvar[j]*res[j] for j in 1:nmult))
  end

  return getobjectivevalue(ms)
end

# Lagrangean Heuristics
function fixbinaries(graph::PlasmoGraph,cat=[:Bin])
  mflat = graph.attributes[:mflat]
  for j in 1:mflat.numCols
    if mflat.colCat[j] in cat
      mflat.colUpper[j] = mflat.colVal[j]
      mflat.colLower[j] = mflat.colVal[j]
    end
  end
  status = solve(mflat)
  if status == :Optimal
    return getobjectivevalue(mflat)
  else
    error("Heuristic model not infeasible or unbounded")
  end
end

function fixintegers(graph::PlasmoGraph)
  fixbinaries(graph,[:Bin,:Int])
end

# Main Function
function  lagrangesolveold(graph::PlasmoGraph;
  update_method=:subgradient,
  max_iterations=100,
  ϵ=0.001,
  α=2,
  UB=5e5,
  LB=-1e5,
  δ=0.8,
  ξ1=0.1,
  ξ2=0,
  λinit=:relaxation,
  solveheuristic=fixbinaries,
  timelimit=360000)

  ########## 0. Initialize ########
  # Start clock
  tic()
  starttime = time()

  # Results outputs
  df = DataFrame(Iter=[],Time=[],α=[],step=[],UB=[],LB=[],Hk=[],Zk=[],Gap=[])
  res = Dict()

  # Get Linkings
  links = getlinkconstraints(graph)
  nmult = length(links)

  # Generate subproblem array
  SP = [getmodel(graph.nodes[i]) for i in 1:length(graph.nodes)]
  for sp in SP
    JuMP.setsolver(sp,graph.solver)
  end
  # Capture objectives
  SPObjectives = [getmodel(graph.nodes[i]).obj for i in 1:length(graph.nodes)]
  sense = SP[1].objSense

  # Generate model for heuristic
  mflat = create_flat_graph_model(graph)
  mflat.solver = graph.solver
  # Restore mflat sense
  if sense == :Max
    mflat.objSense = :Max
    mflat.obj = -mflat.obj
  end
  # Solve realaxation
  solve(mflat,relaxation=true)
  bestbound = getobjectivevalue(mflat)
  if sense == :Max
    UB = bestbound
  else
    LB = bestbound
  end
  debug("Solved LP relaxation with value $bestbound")
  # Set starting λ to the duals of the LP relaxation
  # TODO handle NLP relaxation
  λk = λinit == :relaxation ? mflat.linconstrDuals[end-nmult+1:end] : λk = [0.0 for j in 1:nmult]
  λprev = λk

  # Variables
  θ = 0
  Kprev = [0 for j in 1:nmult]
  i = 0
  direction = nothing
  Zprev = sense == :Max ? UB : LB

  # Master Model (generate only for  planes or bundle methods)
  if update_method in [:cuttingplanes,:bundle]
    ms = Model(solver=graph.solver)
    @variable(ms, η)
    @variable(ms, λ[1:nmult])
    mssense = :Min
    if sense == :Min
      mssense = :Max
    end
    @objective(ms, mssense, η)
  end

  ## <-- Begin Iterations --> ##

  ########## 1. Solve Subproblems ########
  for iter in 1:max_iterations
    debug("*********************")
    debug("*** ITERATION $iter  ***")
    debug("*********************")

    Zk = 0
    improved = false


    # Restore initial objective
    for (j,sp) in enumerate(SP)
      sp.obj = SPObjectives[j]
    end
    # add dualized part
    for l in 1:nmult
      for j in 1:length(links[l].terms.vars)
        var = links[l].terms.vars[j]
        coeff = links[l].terms.coeffs[j]
        var.m.obj += λk[l]*coeff*var
      end
    end

    # Solve
    SP_result = pmap(psolve,SP)
    # Put values back in the graph
    nodedict = getnodes(graph)
    for spd in SP_result
      Zk += spd[:objective]
      getmodel(nodedict[spd[:nodeindex]]).colVal = spd[:values]
    end
    debug("Zk = $Zk")


    ########## 2. Solve Lagrangean Heuristic ########
    mflat.colVal = vcat([getmodel(n).colVal for n in values(nodedict)]...)
    Hk = solveheuristic(mflat)
    debug("Hk = $Hk")


    ########## 3. Check for Bounds Convergence ########
    # Update Bounds
    UBprev = UB
    LBprev = LB
    bestbound_prev = bestbound
    UB = sense == :Max ? min(Zk,UB) : min(Hk,UB)
    LB = sense == :Max ? max(Hk,LB) : max(Zk,LB)

    # Update objective value and calculate gap
    objective = sense == :Max ? LB : UB
    bestbound = sense == :Max ? UB : LB
    graph.objVal = objective
    gap = (UB - LB)/objective

    # Check
    if gap < ϵ
      debug("Converged on bounds to $objective")
      break
    end

    # Increase or restore bestbound improvement counter
    i += bestbound == bestbound_prev ? 1 : -i
    debug("i = $i")

    ########## 3. Check for improvement and update λ ########
    ξ = min(0.1, ξ1 + ξ2*gap)
    improved = sense == :Max ? Zk < Zprev*(1+ξ) : Zk > Zprev*(1+ ξ)
    debug("Compared Zkprev + $(round(ξ*100,1))% = $(Zprev*(1+ξ)) with Zk = $Zk and improved is $improved")
    # Force first step
    if iter == 1
      improved = true
      Zprev = Zk
    end

    # Line search
    # If improvement take step, else reduce α
    if improved
      Zk < Zprev && debug("IMPROVED bound")
      Zprev = Zk
      direction = [getvalue(links[j].terms) for j in 1:nmult]
      λprev = λk
      debug("STEP taken")
    else
      α *= 0.5
    end

    # Restore α
#    if iter % 100 == 0
#       α = 2
#    end
    # Shrink α if stuck
    if iter > 10 && i > 4
      α *= δ
      i = 0
      debug("STUCK, shrink α")
    end

    # Check convergence on α and direction
    if α < 1e-12
      debug("Converged on α = $α")
      break
    end

    normdirection = norm(direction)
    if norm(direction) == 0
      debug("Converged to feasible point")
      break
    end

    # Subgradient update
    difference = sense == :Max ? Zk - LB : UB - Zk

    # Direction correction method method
    μ = direction + θ*Kprev
    if update_method == :subgradient_correction
      if  dot(direction,Kprev) < 0
          θ = normdirection/norm(Kprev)
      else
        θ = 0
      end
    end
    # If the update method is without direction correction Θ = 0 and μ defaults to direction
    step = α*difference/dot(direction,μ)
    λk = λprev - step*μ

    # Check step convergence
    if step < 1e-20
      debug("Converged on step = $step")
      break
    end


    # Update multiplier bounds (Bundle method)
    if update_method == :bundle
      for j in 1:nmult
        setupperbound(λ[j], λprev[j] + step*abs(direction[j]))
        setlowerbound(λ[j], λprev[j] - step*abs(direction[j]))
      end
    end
    # Cutting planes or Bundle
    if update_method in (:cuttingplanes,:bundle)
      if sense == :Max
        @constraint(ms, η >= Zk + sum(λ[j]*direction[j] for j in 1:nmult))
      else
        @constraint(ms, η <= Zk + sum(λ[j]*direction[j] for j in 1:nmult))
      end
      debug("Last cut = $(ms.linconstr[end])")
      if iter > 10
        solve(ms)
        λk = getvalue(λ)
      end
    end

    # Report
    debug("Step = $step")
    debug("α = $α")
    debug("UB = $UB")
    debug("LB = $LB")
    debug("gap = $gap")
    elapsed = round(time()-starttime)
    push!(df,[iter,elapsed,α,step,UB,LB,Hk,Zk,gap])

    res[:Iterations] = iter
    res[:Gap] = gap

    if elapsed > timelimit
       debug("Time Limit exceeded, $elapsed seconds")
       break
    end

  end # Iterations

  # Report
  res[:Objective] = sense == :Min ? UB : LB
  res[:BestBound] = sense == :Min ? LB : UB
  res[:Time] = toc()
  return res, df

end # function

# Parallel model solve function, returns an array of objective values with dimension equal to of elements in the collection for which pmap was applied
function psolve(m::JuMP.Model)
  solve(m)
  d = Dict()
  d[:objective] = getobjectivevalue(m)
  d[:values] = m.colVal
  node = getnode(m)
  for v in values(node.index)
    d[:nodeindex] = v
  end
  #println("Solved node $(d[:nodeindex]) on $(gethostname())")
  return d
end

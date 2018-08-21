"""
  lagrangesolve(graph)
  solve graph using lagrange decomposition
"""
function lagrangesolve(graph;
  max_iterations=10,
  update_method=:subgradient, # :intersectionstep, :ADMM, :cuttingplanes, :bundle
  ϵ=0.001, # ϵ-convergence tolerance
  timelimit=3600,
  α=2, # default subgradient step
  lagrangeheuristic=fixbinaries, # function to calculate the upper bound
  initialmultipliers=:zero, # :relaxation for LP relaxation
  δ = 0.5, # Factor to shrink step when subgradient stuck
  maxnoimprove = 3,
  cpbound=1e6) # Amount of iterations that no improvement is allowed before shrinking step

  ### INITIALIZATION ###
  lgprepare(graph,δ,maxnoimprove,cpbound)
  n = graph.attributes[:normalized]

  if initialmultipliers == :relaxation
    initialrelaxation(graph)
  end
 
  starttime = time()
  λ = graph.attributes[:λ][end]
  res = graph.attributes[:res][end]
  nmult = graph.attributes[:numlinks]
  nodes = [node for node in values(getnodes(graph))]
  if !haskey(graph.attributes, :α)
    graph.attributes[:α] = [1.0α]
  end
  iterval = 0

  ### ITERATIONS ###
  #iterationheader()
  for iter in 1:max_iterations
    variant = iter == 1 ? :default : update_method # Use default version in the first iteration

    iterstart = time()
    # Solve subproblems
    Zk = 0
    x = zeros(graph.attributes[:numlinks],2)
    for node in nodes
       (x,Zkn) = solvenode(node,λ,x,variant)
       Zk += Zkn
    end
    if Zk > graph.attributes[:LB]
      graph.attributes[:LB] = Zk
      graph.attributes[:best_λ] = λ
    end 


    #Zk *= n
    graph.attributes[:steptaken] = true
    # Save info
    push!(graph.attributes[:Zk],Zk)
    push!(graph.attributes[:x],x)
    α = graph.attributes[:α][end]
    push!(graph.attributes[:α], α)

    # If no improvement, increase counter
    if length(graph.attributes[:Zk]) > 1
      if Zk < graph.attributes[:Zk][end-1]
        graph.attributes[:noimprove] += 1
      end
    end 
    # If too many iterations without improvement, decrease :α
    if graph.attributes[:noimprove] >= graph.attributes[:maxnoimprove]
      graph.attributes[:noimprove] = 0
      graph.attributes[:α][end] *= graph.attributes[:δ]
    end

    # TODO Check the order of updates and saving


    # Update residuals
    res = x[:,1] - x[:,2]
    push!(graph.attributes[:res],res)



    # Update multipliers
    println("α = $α")
    (λ, iterval) = updatemultipliers(graph,λ,res,update_method,lagrangeheuristic)
    push!(graph.attributes[:λ], λ)
    if iterval < graph.attributes[:UB]
      graph.attributes[:UB] = iterval
    end 
    push!(graph.attributes[:UB_record], iterval)
    itertime = time() - iterstart
    tstamp = time() - starttime
    saveiteration(graph.attributes[:solution],tstamp,[n*iterval,n*Zk,itertime,tstamp],n)
    printiterationsummary(graph.attributes[:solution],singleline=false)

    # Check convergence
    if graph.attributes[:UB] - graph.attributes[:LB] < ϵ
      graph.attributes[:solution].termination = "Optimal"
      return graph.attributes[:solution]
    end

    # Check time limit
    if tstamp > timelimit
      graph.attributes[:solution].termination = "Time Limit"
      return graph.attributes[:solution]
    end

    #check the value of α
    if graph.attributes[:α][end] < 1e-4
      graph.attributes[:solution].termination = "α too small"
      return graph.attributes[:solution]
    end 

  end
  graph.attributes[:solution].termination = "Max Iterations"
  return graph.attributes[:solution]
end

# Preprocess function
"""
  lgprepare(graph::PlasmoGraph)
  Prepares the graph to apply lagrange decomposition algorithm
"""
function lgprepare(graph::PlasmoGraph, δ=0.5, maxnoimprove=3,cpbound=nothing)
  if haskey(graph.attributes,:preprocessed)
    return true
  end
  n = normalizegraph(graph)
  links = getlinkconstraints(graph)
  nmult = length(links) # Number of multipliers
  graph.attributes[:solution] = Solution(method=:dual_decomposition)
  graph.attributes[:numlinks] = nmult
  graph.attributes[:numx] = nmult / (length(getnodes(graph)) -1)
  graph.attributes[:numnodes] = length(getnodes(graph))
  graph.attributes[:λ] = [zeros(nmult)] # Array{Float64}(nmult)
  graph.attributes[:x] = [] # Linking variables values
  graph.attributes[:res] = [zeros(nmult)] # Residuals
  graph.attributes[:Zk] = [] # Bounds
  #graph.attributes[:mflat] = #create_flat_graph_model(graph)
  #graph.attributes[:mflat].solver = graph.solver
  graph.attributes[:cuts] = []
  graph.attributes[:δ] = δ
  graph.attributes[:noimprove] = 0
  graph.attributes[:maxnoimprove] = maxnoimprove
  graph.attributes[:explore] = []
  graph.attributes[:steptaken] = false
  graph.attributes[:LB_record] = []
  graph.attributes[:UB_record] =[]
  graph.attributes[:LB] = -Inf 
  graph.attributes[:UB] = +Inf 

  # Create Lagrange Master
  ms = Model(solver=graph.solver)
  @variable(ms, η, upperbound=cpbound)
  @variable(ms, λ[1:nmult])
  @objective(ms, Max, η)

  graph.attributes[:lgmaster] = ms

  #map scenario to node & save lagrangean cuts 
  if haskey(getnodes(graph)[1].attributes, :scenario)
    graph.attributes[:node_scenariomap] = Dict()
    graph.attributes[:scenario_nodemap] = Dict()
    for i in keys(getnodes(graph))
      node = getnodes(graph)[i]
      graph.attributes[:node_scenariomap][i] = node.attributes[:scenario]
      graph.attributes[:scenario_nodemap][node.attributes[:scenario]] = i 
      node.attributes[:Zsl] = []
      node.attributes[:μ] = []
    end 
  end 



  # Each node save its initial objective and set a solver if they don't have one
  for n in values(getnodes(graph))
    mn = getmodel(n)
    if mn.solver == JuMP.UnsetSolver()
      mn.solver = graph.solver
    end
    mn.ext[:preobj] = mn.obj
    mn.ext[:multmap] = Dict()
    mn.ext[:varmap] = Dict()
    #map variable name to variable 
    n.attributes[:varname_to_var] = Dict()
    n.attributes[:varname_to_varcat] = Dict()
    for i in 1:length(mn.colNames)
      varname = mn.colNames[i]
      n.attributes[:varname_to_var][varname] = Variable(mn, i)
      n.attributes[:varname_to_varcat][varname] = mn.colCat[i]
    end
  end

  # Maps
  # Multiplier map to know which component of λ to take
  # Varmap knows what values to post where
  for (i,lc) in enumerate(links)
    for j in 1:length(lc.terms.vars)
      var = lc.terms.vars[j]
      var.m.ext[:multmap][i] = (lc.terms.coeffs[j],lc.terms.vars[j], j)
      var.m.ext[:varmap][var] = (i,j)
    end
  end
  for n in values(getnodes(graph))
    mn = getmodel(n)
    println(mn.ext[:multmap])
    println(mn.ext[:varmap])
  end
  graph.attributes[:preprocessed] = true
end

# Solve a single subproblem
function solvenode(node,λ,x,variant=:default)
  m = getmodel(node)
  m.obj = m.ext[:preobj]
  m.ext[:lgobj] = m.ext[:preobj]
  new_μ = Dict()
  println(m.ext[:multmap])
  # Add dualized part to objective function
  for k in keys(m.ext[:multmap])
    coef = m.ext[:multmap][k][1]
    var = m.ext[:multmap][k][2]
    var_name = m.colNames[var.col]
    if !haskey(new_μ, var_name)
      new_μ[var_name] = 0.0
    end 
    new_μ[var_name] += λ[k]*coef
    m.ext[:lgobj] += λ[k]*coef*var
    m.obj += λ[k]*coef*var
    if variant == :ADMM
      j = 3 - m.ext[:varmap][var][2]
      m.obj += 1/2*(coef*var - coef*x[k,j])^2
    end
  end

  push!(node.attributes[:μ], new_μ)

  # Optional: If my residuals are zero, do nothing

  solve(m)

  for l in keys(m.ext[:multmap])
    j = m.ext[:multmap][l][3]
    x[l, j] = getvalue(m.ext[:multmap][l][2])
  end



  objval=MathProgBase.getobjbound(m.internalModel)
  node.attributes[:objective] = objval
  push!(node.attributes[:Zsl], objval)
  node.attributes[:solvetime] = getsolvetime(m)

  return x, objval
end

# Multiplier Initialization
function initialrelaxation(graph)
  if !haskey(graph.attributes,:mflat)
    graph.attributes[:mflat] = create_flat_graph_model(graph)
    graph.attributes[:mflat].solver = graph.solver
  end
  n = graph.attributes[:normalized]
  nmult = graph.attributes[:numlinks]
  mf = graph.attributes[:mflat]
  solve(mf,relaxation=true)
  graph.attributes[:λ][end] = n*mf.linconstrDuals[end-nmult+1:end]
  return getobjectivevalue(mf)
end

function updatemultipliers(graph,λ,res,method,lagrangeheuristic=nothing)
  if method == :subgradient
    subgradient(graph,λ,res,lagrangeheuristic)
  elseif method == :intersectionstep
    intersectionstep(graph,λ,res,lagrangeheuristic)
  elseif method == :probingsubgradient
    probingsubgradient(graph,λ,res,lagrangeheuristic)
  elseif method == :marchingstep
    marchingstep(graph,λ,res,lagrangeheuristic)
  elseif method == :ADMM
    ADMM(graph,λ,res,lagrangeheuristic)
  elseif method == :cuttingplanes
    cuttingplanes(graph,λ,res)
  elseif method == :bundle
    bundle(graph,λ,res,lagrangeheuristic)
  elseif  method == :interactive
    interactive(graph,λ,res,lagrangeheuristic)
  end
end

# Update functions
function subgradient(graph,λ,res,lagrangeheuristic)
  α = graph.attributes[:α][end]
  n = graph.attributes[:normalized]
  bound = n*lagrangeheuristic(graph)
  if bound < graph.attributes[:UB]
    graph.attributes[:UB] = bound
  end
  step = α*abs(graph.attributes[:UB]-graph.attributes[:LB])/(norm(res)^2)
  λ += step*res
  return λ,bound
end

function αeval(αv,graph,bound)
  xv = deepcopy(graph.attributes[:x][end])
  res = graph.attributes[:res][end]
  Zk = graph.attributes[:Zk][end]
  λ = graph.attributes[:λ][end]
  nodes = [node for node in values(getnodes(graph))]
  step = abs(Zk-bound)/(norm(res)^2)
  zk = 0
  for node in nodes
     (xv,Zkn) = solvenode(node,λ+αv*step*res,xv,:default)
     zk += Zkn
  end
  return zk
end

function αexplore(graph,bound)
  df = graph.attributes[:explore]
  n = graph.attributes[:normalized]
  z = Float64[]
  for α in 0:0.1:2
    push!(z,n*αeval(α,graph,bound))
  end
  push!(df,z)
end

function probingsubgradient(graph,λ,res,lagrangeheuristic,α=graph.attributes[:α][end],Δ=0.01;exhaustive=false)
  res = graph.attributes[:res][end]
  Zk = graph.attributes[:Zk][end]
  n = graph.attributes[:normalized]
  bound = n*lagrangeheuristic(graph)
  step = abs(Zk-bound)/(norm(res)^2)
  # First point
  α1 = 0
  z1 = Zk

  # Second point
  α2 = α
  z2 = αeval(α2,graph,bound)

  if (z2 - z1)/abs(z1) >= Δ
    # If second point gives an increase of more than Δ%, take it.
    return (λ += α2*step*res), bound
  elseif abs((z2 - z1)/z1) < Δ
    # If second point is similar to the first one, take the midpoint
    return (λ += α2/2*step*res), bound
  else
    # If second point is larger than first, take quarter-step
    if exhaustive
      return probingsubgradient(graph,λ,res,lagrangeheuristic,α/4,Δ;exhaustive=true)
    else
      return (λ += α2/4*step*res), bound
    end
  end
end

function marchingstep(graph,λ,res,lagrangeheuristic,α=graph.attributes[:α][end],Δ=0.1α)
  res = graph.attributes[:res][end]
  Zk = graph.attributes[:Zk][end]
  n = graph.attributes[:normalized]
  bound = n*lagrangeheuristic(graph)
  step = abs(Zk-bound)/(norm(res)^2)
  # First point
  α1 = 0
  z1 = Zk
  zs_1 = z1
  αs_1 = α1

  for αs in α1+Δ:Δ:α
    zs = αeval(αs,graph,bound)
    if zs < zs_1
      return (λ += αs_1*step*res), bound
    end
    zs_1 = zs
    αs_1 = αs
  end

  return (λ += α*step*res), bound
end

@require Gaston begin
function interactive(graph,λ,res,lagrangeheuristic)
  α = graph.attributes[:α][end]
  n = graph.attributes[:normalized]
  bound = n*lagrangeheuristic(graph)
  Zk = graph.attributes[:Zk][end]
  αexplore(graph,bound)
  plot(0:0.1:2,graph.attributes[:explore][end])
  print("α = ")
  α = parse(Float64,readline(STDIN))
  step = α*abs(Zk-bound)/(norm(res)^2)
  λ += step*res
  return λ,bound
end
end

function intersectionstep(graph,λ,res,lagrangeheuristic,α=graph.attributes[:α][end],Δ=0.01,ϵ=0.001)
  res = graph.attributes[:res][end]
  Zk = graph.attributes[:Zk][end]
  n = graph.attributes[:normalized]
  bound = n*lagrangeheuristic(graph)
  step = abs(Zk-bound)/(norm(res)^2)
  # First curve
  αa0 = 0
  za0 = Zk
  αa1 = Δ
  za1 = αeval(αa1,graph,bound)
  ma = (za1 - za0)/(αa1 - αa0)
  if abs(ma) < ϵ
    return λ,bound
  elseif ma < 0
    warn("First slope decreasing. This might be an indication that there is no improvement in the direction chosen")
    return λ,bound
  end

  # Second curve
  αb0 = α
  zb0 = αeval(αb0,graph,bound)
  αb1 = αb0 - Δ
  zb1 = αeval(αb1,graph,bound)
  mb = (zb1 - zb0)/(αb1 - αb0)
  println("ma = $ma")
  println("mb = $mb")
  if abs(mb) < ϵ
    return λ + α*step*res, bound
  end
  # Check different Sign
  if sign(ma)>0 && sign(mb)>0
    return λ + α*step*res, bound
  end
  # Find intersection
  αinter = (za0 - zb0 + αb0*mb)/(mb - ma)
  λ += αinter*step*res
  return λ,bound
end

function cuttingplanes(graph,λ,res)
  ms = graph.attributes[:lgmaster]
  Zk = graph.attributes[:Zk][end]
  nmult = graph.attributes[:numlinks]

  λvar = getindex(ms, :λ)
  η = getindex(ms,:η)

  cut = @constraint(ms, η <= Zk + sum(λvar[j]*res[j] for j in 1:nmult))
  push!(graph.attributes[:cuts], cut)

  solve(ms)
  return getvalue(λvar), getobjectivevalue(ms)
end

function bundle(graph,λ,res,lagrangeheuristic)
  α = graph.attributes[:α][end]
  bound = lagrangeheuristic(graph)
  Zk = graph.attributes[:Zk][end]
  ms = graph.attributes[:lgmaster]
  λvar = getindex(ms, :λ)
  step = α*abs(Zk-bound)/(norm(res)^2)
  setlowerbound.(λvar,λ-step*abs.(res))
  setupperbound.(λvar,λ+step*abs.(res))

  cuttingplanes(graph,λ,res)
end

function ADMM(graph,λ,res,lagrangeheuristic)
  bound = lagrangeheuristic(graph)
  λ += res/norm(res)
  return λ,bound
end

# Lagrangean Heuristics
function fixbinaries(graph::PlasmoGraph,cat=[:Bin])
  if !haskey(graph.attributes,:mflat)
    graph.attributes[:mflat] = create_flat_graph_model(graph)
  end
  n = graph.attributes[:normalized]
  mflat = graph.attributes[:mflat]
  mflat.solver = graph.solver
  mflat.colVal = vcat([getmodel(n).colVal for n in values(getnodes(graph))]...)
  for j in 1:mflat.numCols
    if mflat.colCat[j] in cat
      mflat.colUpper[j] = mflat.colVal[j]
      mflat.colLower[j] = mflat.colVal[j]
    end
  end
  status = solve(mflat)
  if status == :Optimal
    return n*getobjectivevalue(mflat)
  else
    error("Heuristic model not infeasible or unbounded")
  end
end

function fixintegers(graph::PlasmoGraph)
  fixbinaries(graph,[:Bin,:Int])
end

function nearest_scenario(graph::PlasmoGraph)
  avg_x = Dict()
  x_all_scenarios = []
  for i in 1:graph.attributes[:numnodes]
    push!(x_all_scenarios, Dict())
  end 
  orig_ub = Dict()
  orig_lb = Dict()
  #calculate the average of x over all scenarios
  j=1
  for node in values(getnodes(graph))
    m  = getmodel(node)
    # println(keys(m.ext[:varmap]))
    for var in keys(m.ext[:varmap])
      var_name = m.colNames[var.col]
      avg_x[var_name] = 0.0
    end
    # println(keys(m.ext[:varmap]))
    for var in keys(m.ext[:varmap])
      var_name = m.colNames[var.col]
      avg_x[var_name] += node.attributes[:prob] * getvalue(var)
      x_all_scenarios[j][var_name] = getvalue(var)
      if j == 1
        orig_ub[var_name] = getupperbound(var)
        orig_lb[var_name] = getlowerbound(var)
      end
    end
    j += 1
  end
  println("x_all_scenarios")
  println(x_all_scenarios)
  #find the nearest scenario 
  nearest_x = Dict()
  min_distance = +Inf  
  for j in 1:length(graph.attributes[:numnodes])
    temp_distance = 0 
    for var_name in keys(avg_x)
      temp_distance += ((avg_x[var_name] - x_all_scenarios[j][var_name]) / (orig_ub[var_name]  - orig_lb[var_name]))^2
    end
    if temp_distance < min_distance
      min_distance = temp_distance
      for var_name in keys(avg_x)
        nearest_x[var_name] = x_all_scenarios[j][var_name]
      end
    end 
  end


  #fix x to nearest_x and re-solve all subproblems
  Zk = 0 
  for node in values(getnodes(graph))
    m = getmodel(node)
    for var in keys(m.ext[:varmap])
      var_name = var.m.colNames[var.col]
      setupperbound(var, nearest_x[var_name])
      setlowerbound(var, nearest_x[var_name])
    end 
    status = solve(m)
    # if status != :Optimal
    #   error("upper bound subproblem not solved to optimality")
    # end 
    Zk += getobjectivevalue(m)
    #restore bounds after solve 
    for var in keys(m.ext[:varmap])
      var_name = var.m.colNames[var.col]
      setupperbound(var, orig_ub[var_name])
      setlowerbound(var, orig_lb[var_name])
    end
  end
  if Zk < graph.attributes[:UB]
    graph.attributes[:best_feasible_x]  = nearest_x
  end 
  return Zk

end


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

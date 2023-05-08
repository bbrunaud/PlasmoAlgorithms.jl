"""
  lagrangesolve(graph)
  solve graph using lagrange decomposition
"""
function lagrangeoptimize!(graph;
  max_iterations=10,
  update_method=:subgradient, #probingsubgradient
  ϵ=0.001, # ϵ-convergence tolerance
  timelimit=3600,
  α=2, # default subgradient step
  lagrangeheuristic=fixbinaries, # function to calculate the upper bound
  initialmultipliers=:zero, # :relaxation for LP relaxation
  δ = 0.5, # Factor to shrink step when subgradient stuck
  maxnoimprove = 3,
  cpbound=1e6,
  verbose=true,
  singleline = true) # Amount of iterations that no improvement is allowed before shrinking step

  ### INITIALIZATION ###
  lgprepare(graph, δ=δ, maxnoimprove=maxnoimprove,cpbound=cpbound)
  n = getattribute(graph,:normalized)

  if initialmultipliers == :relaxation
    initialrelaxation(graph)
  end

  starttime = time()
  s = Solution(method=:dual_decomposition)
  λ = getattribute(graph,:λ)[end]
  x = getattribute(graph,:x)[end]
  res = getattribute(graph,:res)[end]
  nmult = getattribute(graph,:numlinks)
  nodes = [node for node in values(getnodes(graph))]
  setattribute(graph,:α, [1.0α])
  iterval = 0

  ### ITERATIONS ###
  for iter in 1:max_iterations
    variant = iter == 1 ? :default : update_method # Use default version in the first iteration

    iterstart = time()
    # Solve subproblems
    Zk = 0
    for node in nodes
       (x,Zkn) = solvenode(graph,node,λ,x,variant)
       Zk += Zkn
    end
    setattribute(graph, :steptaken, true)

    # If no improvement, increase counter
    if Zk < getattribute(graph, :Zk)[end]
      setattribute(graph, :noimprove, getattribute(graph, :noimprove) + 1)
    end
    # If too many iterations without improvement, decrease :α
    if getattribute(graph, :noimprove) >= getattribute(graph, :maxnoimprove)
      setattribute(graph, :noimprove, 0)
      getattribute(graph, :α)[end] *= getattribute(graph, :δ)
    end
    # Save info
    push!(getattribute(graph, :Zk),Zk)
    push!(getattribute(graph, :x),x)
    α = getattribute(graph, :α)[end]
    push!(getattribute(graph, :α), α)

    # Update residuals
    res = x[:,1] - x[:,2]
    push!(getattribute(graph, :res),res)

    # Save iteration data
    itertime = time() - iterstart
    tstamp = time() - starttime
    saveiteration(s,tstamp,[n*iterval,n*Zk,itertime,tstamp],n)
    if verbose 
      iter == 1 && singleline && printheader(graph, :Lagrange)
      printiterationsummary(s,singleline=singleline)
    end
    
    # Check convergence
    if norm(res) < ϵ
      s.termination = "Optimal"
      return s
    end

    # Check time limit
    if tstamp > timelimit
      s.termination = "Time Limit"
      return s
    end

    # Update multipliers
    (λ, iterval) = updatemultipliers(graph,λ,res,update_method,lagrangeheuristic)
    push!(getattribute(graph, :λ), λ)
    # Update iteration time
    s.itertime[end] = time() - iterstart
  end
  s.termination = "Max Iterations"
  return s
end

# Preprocess function
"""
  lgprepare(graph::PlasmoGraph)
  Prepares the graph to apply lagrange decomposition algorithm
"""
function lgprepare(graph::OptiGraph; δ=0.5, maxnoimprove=3,cpbound=1e6)
  if hasattribute(graph,:preprocessed)
    return true
  end
  n = normalizegraph(graph)
  links = getlinkconstraints(graph)
  nmult = length(links) # Number of multipliers
  setattribute(graph, :numlinks, nmult)
  setattribute(graph, :λ, [zeros(nmult)]) # Array{Float64}(nmult)
  setattribute(graph, :x, [zeros(nmult,2)]) # Linking variables values
  setattribute(graph, :res, [zeros(nmult)]) # Residuals
  setattribute(graph, :Zk, [0.0]) # Bounds
  setattribute(graph, :cuts, [])
  setattribute(graph, :δ, δ)
  setattribute(graph, :noimprove, 0)
  setattribute(graph, :maxnoimprove, maxnoimprove)
  setattribute(graph, :explore, [])
  setattribute(graph, :steptaken, false)
  numnodes = length(getnodes(graph))
  setattribute(graph, :numnodes, numnodes)
  setattribute(graph, :cutdata, Dict(i => CutData[] for i in 1:numnodes))

  # Create Lagrange Master
  ms = isnothing(graph.optimizer) ? Model() : Model(graph.optimizer)
  #set_silent(ms)
  @variable(ms, η[1:numnodes], upper_bound=cpbound)
  @variable(ms, λ[1:nmult])
  @objective(ms, Max, sum(η))

  setattribute(graph, :lgmaster, ms)

  # Each node save its initial objective and set a solver if they don't have one
  for n in values(getnodes(graph))
    mn = getmodel(n)
    mn.ext[:preobj] = objective_function(mn)
    mn.ext[:multmap] = Dict()
    mn.ext[:varmap] = Dict()
  end

  # Maps
  # Multiplier map to know which component of λ to take
  # Varmap knows what values to post where
  for (i,lc) in enumerate(links)
    for (j,term) in enumerate(lc.func.terms)
      var = term[1]
      coeff = term[2]
      owner_model(var).ext[:multmap][i] = (coeff, var)
      owner_model(var).ext[:varmap][var] = (i,j)
    end
  end

  setattribute(graph, :preprocessed, true)
end

# Solve a single subproblem
function solvenode(graph, node,λ,x,variant=:default)
  m = getmodel(node)
  # Restore objective function
  set_objective_function(m, m.ext[:preobj])  
  m.ext[:lgobj] = m.ext[:preobj]
  # Add dualized part to objective function
  for k in keys(m.ext[:multmap])
    coef = m.ext[:multmap][k][1]
    var = m.ext[:multmap][k][2]
    m.ext[:lgobj] += λ[k]*coef*var
    set_objective_function(m, objective_function(m) + λ[k]*coef*var)
    if variant == :ADMM
      j = 3 - m.ext[:varmap][var][2]
      set_objective_function(m, objective_function(m) + 1/2*(coef*var - coef*x[k,j])^2)
    end
  end
  # Solve
  optimize!(m)
  # Pass output
  for v in keys(m.ext[:varmap])
    val = value(v)
    x[m.ext[:varmap][v]...] = val
  end
  objval = value(m.ext[:lgobj])
  setattribute(node, :objective, objval)
  node[:solvetime] =  solve_time(m)

  # Push Cut Data
  nodeindex = graph[node]
  preobjval = value(m.ext[:preobj])
  λcomponent = Int64[]
  coeffs = Float64[]
  xk = Float64[]
  for k in keys(m.ext[:multmap])
    coeff = m.ext[:multmap][k][1]
    var = m.ext[:multmap][k][2]
    push!(λcomponent, k)
    push!(coeffs, coeff)
    push!(xk, value(var))
  end
  cutdataarray = getattribute(graph,:cutdata)
  push!(cutdataarray[nodeindex], LagrangeCutData(preobjval, λcomponent, coeffs, xk))

  return x, objval
end

# Multiplier Initialization
function initialrelaxation(graph)
  if !hasattribute(graph,:mflat)
    optinode, = combine(graph)
    optimodel = getmodel(optinode)
    JuMP.set_optimizer(optimodel, graph.optimizer)
    setattribute(graph, :mflat, optimodel)
  end
  n = getattribute(graph , :normalized)
  nmult = getattribute(graph , :numlinks)
  mf = getattribute(graph , :mflat)
  unrelax = relax_integrality(mf)
#  set_silent(mf)
  optimize!(mf)
  λ = getattribute(graph , :λ)
  equalities = all_constraints(mf, AffExpr, MOI.EqualTo{Float64})
  links = equalities[1:nmult]
  λ[end] = n*dual.(links)
#  unrelax()
  return objective_value(mf)
end


function updatemultipliers(graph,λ,res,method,lagrangeheuristic=nothing)
  if method == :subgradient
    subgradient(graph,λ,res,lagrangeheuristic)
  elseif method == :probingsubgradient
    probingsubgradient(graph,λ,res,lagrangeheuristic)
  elseif method == :cuttingplanes
    cuttingplanes(graph)
  elseif  method == :interactive
    interactive(graph,λ,res,lagrangeheuristic)
  end
end

# Update functions
function subgradient(graph,λ,res,lagrangeheuristic)
  α = getattribute(graph , :α)[end]
  n = getattribute(graph , :normalized)
  bound = n*lagrangeheuristic(graph)
  Zk = getattribute(graph , :Zk)[end]
  step = α*abs(Zk-bound)/(norm(res)^2)
  λ += step*res
  return λ,bound
end

function αeval(αv,graph,bound)
  xv = deepcopy(getattribute(graph , :x)[end])
  res = getattribute(graph , :res)[end]
  Zk = getattribute(graph , :Zk)[end]
  λ = getattribute(graph , :λ)[end]
  nodes = [node for node in getnodes(graph)]
  step = abs(Zk-bound)/(norm(res)^2)
  zk = 0
  for node in nodes
     (xv,Zkn) = solvenode(graph, node,λ+αv*step*res,xv,:default)
     zk += Zkn
  end
  return zk
end

function αexplore(graph,bound)
  df = getattribute(graph , :explore)
  n = getattribute(graph , :normalized)
  z = Float64[]
  for α in 0:0.1:2
    push!(z,n*αeval(α,graph,bound))
  end
  push!(df,z)
end

function probingsubgradient(graph,λ,res,lagrangeheuristic,α=getattribute(graph , :α)[end],Δ=0.01;exhaustive=false)
  res = getattribute(graph , :res)[end]
  Zk = getattribute(graph , :Zk)[end]
  n = getattribute(graph , :normalized)
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

function interactive(graph,λ,res,lagrangeheuristic)
  α = getattribute(graph , :α)[end]
  n = getattribute(graph , :normalized)
  bound = n*lagrangeheuristic(graph)
  Zk = getattribute(graph , :Zk)[end]
  αexplore(graph,bound)
  plot(0:0.1:2,getattribute(graph , :explore)[end])
  print("α = ")
  α = parse(Float64,readline(STDIN))
  step = α*abs(Zk-bound)/(norm(res)^2)
  λ += step*res
  return λ,bound
end

function cuttingplanes(graph)
  ms = getattribute(graph, :lgmaster)
  nmult = getattribute(graph, :numlinks)

  λvar = ms[:λ]
  η = ms[:η]

  cutdataarray = getattribute(graph,:cutdata)
  for node in getnodes(graph)
    nodeindex = graph[node]
    while length(cutdataarray[nodeindex]) > 0
      cd = pop!(cutdataarray[nodeindex])
      cut = @constraint(ms, η[nodeindex] <= cd.zk + sum(cd.coeffs[j]*λvar[cd.λc[j]]*cd.xk[j] for j in 1:length(cd.λc)))
      push!(getattribute(graph , :cuts), cut)
    end
  end

  optimize!(ms)
  return value.(λvar), objective_value(ms)
end

#=
# Standard Lagrangean Heuristics
function fixbinaries(graph::OptiGraph,cat=[:Bin])
  if !hasattribute(graph,:mflat)
    setattribute(graph, :mflat, create_jump_graph_model(graph))
    getattribute(graph, :mflat).solver = getsolver(graph)
  end
  n = getattribute(graph, :normalized)
  mflat = getattribute(graph , :mflat)
  mflat.colVal = vcat([getmodel(n).colVal for n in getnodes(graph)]...)
  for j in 1:mflat.numCols
    if mflat.colCat[j] in cat
      mflat.colUpper[j] = mflat.colVal[j]
      mflat.colLower[j] = mflat.colVal[j]
    end
  end
  println("Solving Heuristic with Fixing Binary/Integer Variables")
  status = optimize!(mflat)
  if status == :Optimal
    return n*objective_value(mflat)
  else
    error("Heuristic model not infeasible or unbounded")
  end
end

function fixintegers(graph::OptiGraph)
  fixbinaries(graph,[:Bin,:Int])
end
=#

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
  combinationdef=[]) # Amount of iterations that no improvement is allowed before shrinking step

  lgprepare(graph,δ,maxnoimprove)
  n = graph.attributes[:normalized]

  if initialmultipliers == :relaxation
    initialrelaxation(graph)
  end

  starttime = time()
  s = Solution(method=:dual_decomposition)
  λ = graph.attributes[:λ][end]
  x = graph.attributes[:x][end]
  res = graph.attributes[:res][end]
  nmult = graph.attributes[:numlinks]
  nodes = [node for node in values(getnodes(graph))]
  graph.attributes[:α] = [1.0α]
  iterval = 0


  for iter in 1:max_iterations
    variant = iter == 1 ? :default : update_method # Use default version in the first iteration

    iterstart = time()
    # Solve subproblems
    Zk = 0
    for node in nodes
       (x,Zkn) = solvenode(node,λ,x,variant)
       Zk += Zkn
    end
    Zk *= n
    if Zk > graph.attributes[:Zk][end]
      graph.attributes[:noimprove] += 1
    end
    if graph.attributes[:noimprove] >= graph.attributes[:maxnoimprove]
      graph.attributes[:noimprove] = 0
      graph.attributes[:α][end] *= graph.attributes[:δ]
    end
    push!(graph.attributes[:Zk],Zk)
    push!(graph.attributes[:x],x)


    # Update residuals
    res = x[:,1] - x[:,2]
    push!(graph.attributes[:res],res)

    itertime = time() - iterstart
    tstamp = time() - starttime
    saveiteration(s,tstamp,[iterval,Zk,itertime,tstamp],n)
    iterationsummary(s)

    # Check convergence
    if norm(res) < ϵ
      s.termination = "Optimal"
      return s
    end

    # Update multipliers
    push!(graph.attributes[:α], α)
    (λ, iterval) = updatemultipliers(graph,λ,res,update_method,lagrangeheuristic)
    push!(graph.attributes[:λ], λ)
    # Save summary
  end
  s.termination = "Max Iterations"
  return s
end

# Preprocess function
"""
  lgprepare(graph::PlasmoGraph)
  Prepares the graph to apply lagrange decomposition algorithm
"""
function lgprepare(graph::PlasmoGraph, δ=0.5, maxnoimprove=3)
  if haskey(graph.attributes,:preprocessed)
    return true
  end
  n = normalizegraph(graph)
  links = getlinkconstraints(graph)
  nmult = length(links) # Number of multipliers
  graph.attributes[:numlinks] = nmult
  graph.attributes[:λ] = [zeros(nmult)] # Array{Float64}(nmult)
  graph.attributes[:x] = [zeros(nmult,2)] # Linking variables values
  graph.attributes[:res] = [zeros(nmult)] # Residuals
  graph.attributes[:Zk] = [0.0] # Residuals
  graph.attributes[:mflat] = create_flat_graph_model(graph)
  graph.attributes[:mflat].solver = graph.solver
  graph.attributes[:cuts] = []
  graph.attributes[:δ] = δ
  graph.attributes[:noimprove] = 0
  graph.attributes[:maxnoimprove] = maxnoimprove
  graph.attributes[:df] = []

  # Create Lagrange Master
  ms = Model(solver=graph.solver)
  @variable(ms, η, upperbound=1e-6)
  @variable(ms, λ[1:nmult])
  @objective(ms, Max, η)

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
function solvenode(node,λ,x,variant=:default)
  m = getmodel(node)
  m.obj = m.ext[:preobj]
  m.ext[:lgobj] = m.ext[:preobj]
  # Add dualized part to objective function
  for k in keys(m.ext[:multmap])
    coef = m.ext[:multmap][k][1]
    var = m.ext[:multmap][k][2]
    m.ext[:lgobj] += λ[k]*coef*var
    m.obj += λ[k]*coef*var
    if variant == :ADMM
      j = 3 - m.ext[:varmap][var][2]
      m.obj += 1/2*(coef*var - coef*x[k,j])^2
    end
  end

  # Optional: If my residuals are zero, do nothing

  solve(m)
  for v in keys(m.ext[:varmap])
    val = getvalue(v)
    x[m.ext[:varmap][v]...] = val
  end

  objval = getvalue(m.ext[:lgobj])
  node.attributes[:objective] = objval
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
  elseif method == :bettersubgradient
    bettersubgradient(graph,λ,res,lagrangeheuristic)
  elseif method == :fastsubgradient
    fastsubgradient(graph,λ,res,lagrangeheuristic)
  elseif method == :marchingstep
    marchingstep(graph,λ,res,lagrangeheuristic)
  elseif method == :ADMM
    ADMM(graph,λ,res,lagrangeheuristic)
  elseif method == :cuttingplanes
    cuttingplanes(graph,λ,res)
  elseif method == :bundle
    bundle(graph,λ,res,lagrangeheuristic)
  end
end

# Update functions
function subgradient(graph,λ,res,lagrangeheuristic)
  α = graph.attributes[:α][end]
  bound = lagrangeheuristic(graph)
  Zk = graph.attributes[:Zk][end]
  #αexplore(graph,bound)
  step = α*abs(Zk-bound)/(norm(res)^2)
  λ += step*res
  return λ,bound
end

function αeval(αv,graph,bound)
  xv = deepcopy(graph.attributes[:x][end])
  res = graph.attributes[:res][end]
  Zk = graph.attributes[:Zk][end]
  n = graph.attributes[:normalized]
  λ = graph.attributes[:λ][end]
  nodes = [node for node in values(getnodes(graph))]
  step = abs(Zk-bound)/(norm(res)^2)
  zk = 0
  for node in nodes
     (xv,Zkn) = solvenode(node,λ+αv*step*res,xv,:default)
     zk += Zkn
  end
  zk *= n
  return zk
end

function αexplore(graph,bound)
  df = graph.attributes[:df]
  z = []
  for α in -2:0.1:2
    push!(z,αeval(α,graph,bound))
  end
  push!(df,z)
end


function intersectionstep(graph,λ,res,lagrangeheuristic,α=graph.attributes[:α][end],Δ=0.01,ϵ=0.001)
  res = graph.attributes[:res][end]
  Zk = graph.attributes[:Zk][end]
  bound = lagrangeheuristic(graph)
  step = abs(Zk-bound)/(norm(res)^2)
  #αexplore(graph,bound)
  # First curve
  αa0 = 0
  za0 = Zk
  αa1 = Δ
  za1 = αeval(αa1,graph,bound)
  ma = (za1 - za0)/(αa1 - αa0)
  if abs(ma) < ϵ
    return λ,bound
  elseif ma > 0
    #return λ - α*step*res, bound #
    α = -α
    Δ = -Δ
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
  if sign(ma)<0 && sign(mb)<0
    return λ + α*step*res, bound
  end
  # Find intersection
  αinter = (za0 - zb0 + αb0*mb)/(mb - ma)
  λ += αinter*step*res
  return λ,bound
end

function fastsubgradient(graph,λ,res,lagrangeheuristic,α=graph.attributes[:α][end],Δ=0.01)
  res = graph.attributes[:res][end]
  Zk = graph.attributes[:Zk][end]
  bound = lagrangeheuristic(graph)
  step = abs(Zk-bound)/(norm(res)^2)
  # First point
  α1 = 0
  z1 = Zk

  # Second point
  α2 = α
  z2 = αeval(α2,graph,bound)

  if (z1 - z2)/z1 > Δ
    # If second point gives a decrease of more than Δ%, take it.
    return (λ += α2*step*res), bound
  elseif abs((z1 - z2)/z1) < Δ
    # If second point is similar to the first one, take the midpoint
    return (λ += α2/2*step*res), bound
  else
    # If second point is larger than first, take quarter-step
    return (λ += α2/4*step*res), bound
  end
end

function bettersubgradient(graph,λ,res,lagrangeheuristic,α=graph.attributes[:α][end],Δ=0.01)
  res = graph.attributes[:res][end]
  Zk = graph.attributes[:Zk][end]
  bound = lagrangeheuristic(graph)
  step = abs(Zk-bound)/(norm(res)^2)
  # First point
  α1 = 0
  z1 = Zk

  # Second point
  α2 = α
  debug("α = $α")
  z2 = αeval(α2,graph,bound)

  if (z1 - z2)/z1 > Δ
    # If second point gives a decrease of more than Δ%, take it.
    return (λ += α2*step*res), bound
  elseif abs((z1 - z2)/z1) < Δ
    # If second point is similar to the first one, take the midpoint
    return (λ += α2/2*step*res), bound
  end

  # Third point
  α3 = α/4
  z3 = αeval(α3,graph,bound)
  if z3 < z1
    # If the third point gives decrease, take the step.
    return (λ += α3*step*res), bound
  else
    # If the third point does not give decrease, take a midpoint between them
    return (λ += α3/2*step*res), bound
  end
end

function marchingstep(graph,λ,res,lagrangeheuristic,α=graph.attributes[:α][end],Δ=0.1α)
  res = graph.attributes[:res][end]
  Zk = graph.attributes[:Zk][end]
  bound = lagrangeheuristic(graph)
  step = abs(Zk-bound)/(norm(res)^2)
  # First point
  α1 = 0
  z1 = Zk
  zs_1 = z1
  αs_1 = α1

  for αs in α1+Δ:Δ:α
    zs = αeval(αs,graph,bound)
    if zs > zs_1
      return (λ += αs_1*step*res), bound
    end
    zs_1 = zs
    αs_1 = αs
  end

  return (λ += α*step*res), bound
end


function ADMM(graph,λ,res,lagrangeheuristic)
  bound = lagrangeheuristic(graph)
  λ += res/norm(res)
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

# Lagrangean Heuristics
function fixbinaries(graph::PlasmoGraph,cat=[:Bin])
  if !haskey(graph.attributes,:mflat)
    graph.attributes[:mflat] = create_flat_graph_model(graph)
  end
  n = graph.attributes[:normalized]
  mflat = graph.attributes[:mflat]
  mflat.solver = graph.solver
  mflat.colVal = vcat([getmodel(n).colVal for n in values(getnodes(g))]...)
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

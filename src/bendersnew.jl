"""
bendersolve
"""
function bendersolve(graph::Plasmo.PlasmoGraph; max_iterations::Int64=10, cuts::Array{Symbol,1}=[:LP], ϵ=1e-5,UB=Inf,UBupdatefrequency=1)
  starttime = time()
  tmpdir = "/tmp/RootNode" # mktempdir()
  s = Solution(method=:benders)
  updatebound = true

  bdprepare(graph)
  n = graph.attributes[:normalized]

  solve(graph.attributes[:mflat],relaxation=true)
  LB = n*getobjectivevalue(graph.attributes[:mflat])
  UB = n*UB

  for i in 1:max_iterations
    tic()
    updatebound = i-1 % UBupdatefrequency
    LB,UB = forwardstep(graph, cuts, updatebound)

    if abs(UB-LB) < ϵ
      s.termination = "Optimal"
      return s
    end

    # Check time limit
    tstamp = starttime - time()
    if tstamp > timelimit
      s.termination = "Time Limit"
      return s
    end
    saveiteration(s,tstamp,[n*UB,n*LB,toc(),tstamp],n)
    printiterationsummary(s,singleline=false)
  end

  s.termination = "Max Iterations"
  return s
end

function forwardstep(graph::PlasmoGraph, cuts::Array{Symbol,1}, updatebound::Bool)
  levels = getlevels(graph)
  numlevels = length(levels)
  for level in 1:numlevels
    nodeslevel = levels[level]
    for node in nodeslevel
      solvenode(node,graph,cuts,updatebound)
    end
  end
end

function solvenode(node::PlasmoNode, graph::PlasmoGraph, cuts::Array{Symbol,1}, updatebound::Bool)
  nodenumber = getnodenumber(node)
  numinneighbors = length(in_neighbors(graph,node))
  numoutneighbors = length(out_neighbors(graph,node))
  x = getx(graph)
  λ = getλ(graph)
  nodebound = getnodebounds(graph)

  # 1. Add cuts
  if numoutneighbors > 0   # Not a leaf node
    generatebenderscut(node,nodebound,x)
  end
  # 2. Take x
  if numinneighbors > 0    # Not a root node
    takex(node,x)
  end
  # 3. solve
  if :LP in cuts
    xnode, λnode, nodebound = solvelprelaxation(node)
  end
  if :Root in cuts
    xnode, λnode, nodebound = solverootrelaxation(node)
  end
  if updatebound
    xnode, znode = solvenode(node)
  end
  # 4. put x
  if numoutneighbors > 0   # Not a leaf node
    putx(node,xnode,x)
  end
  # 5. put λ and nodebound
  if numinneighbors > 0    # Not a root node
    putλ(node,λnode,λ)
  end
end

function solvelprelaxation(node)
  model = getmodel(node)
  status = solve(model, relaxation = true)

  @assert status == :Optimal

  dualCon = dualMap[node]
  nodeoutvars = outxMap[node]

  λnode = getdual(dualCon)
  nodebound = getobjectivevalue(model)
  xnode = getvalue(nodeoutvars)

  return  xnode, λnode, nodebound
end

function solverootrelaxation(node)
  sp = getmodel(node)
  lpfile = joinpath(tmpdir,"nodemodel.lp")
  writeLP(sp,lpfile)
  run(`cpxgetroot nodemodel.lp 0 1`)
  lp = Model(solver=sp.solver)
  turnoffpresolve(lp)
  lp.internalModel = MathProgBase.LinearQuadraticModel(lp.solver)
  MathProgBase.loadproblem!(lp.internalModel,"node0.lp")
  MathProgBase.optimize!(lp.internalModel)
  run(`mv node0.lp $tmpdir/`)
  dualCon = dualMap[node]
  nodeoutvars = outxMap[node]

  rootduals = MathProgBase.getconstrduals(lp.internalModel)
  sp.linconstrDuals = MathProgBase.getconstrduals(lp.internalModel)[1:length(sp.linconstrDuals)]

  λnode = getdual(dualCon)
  nodebound = MathProgBase.getobjval(lp.internalModel)
  xnode = getvalue(nodeoutvars)

  return  xnode, λnode, nodebound
end

function generatebenderscut(node::PlasmoNode, nodebound::Real, x)
  println("In generation function")
end

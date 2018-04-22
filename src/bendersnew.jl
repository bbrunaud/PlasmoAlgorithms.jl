
abstract type CutData end

struct BendersCutData <: CutData
  θk
  λk
  xk
end

struct LLIntegerCutData <: CutData
  θlb
  yk
end

struct IntegerCutData <: CutData
  yk
end


"""
bendersolve
"""
function bendersolve(graph::Plasmo.PlasmoGraph; max_iterations::Int64=10, cuts::Array{Symbol,1}=[:LP], ϵ=1e-5,UBupdatefrequency=1,timelimit=3600)
  starttime = time()
  global tmpdir = "/tmp/RootNode" # mktempdir()
  s = Solution(method=:benders)
  updatebound = true

  bdprepare(graph)
  n = graph.attributes[:normalized]

  mf = graph.attributes[:mflat]
  solve(mf,relaxation=true)
  LB = getobjectivevalue(graph.attributes[:mflat])
  UB = Inf

  # Set bound to root node
  rootnode = graph.attributes[:roots][1]
  rootmodel = getmodel(rootnode)
  @constraint(rootmodel, rootmodel.obj.aff >= LB)

  # Begin iterations
  for i in 1:max_iterations
    tic()
    updatebound = ((i-1) % UBupdatefrequency) == 0
    LB,UB = forwardstep(graph, cuts, updatebound)

    tstamp = starttime - time()
    if n == 1
      saveiteration(s,tstamp,[UB,LB,toc(),tstamp],n)
    else
      saveiteration(s,tstamp,[n*LB,n*UB,toc(),tstamp],n)
    end
    printiterationsummary(s,singleline=false)

    if abs(UB-LB) < ϵ
      s.termination = "Optimal"
      return s
    end

    # Check time limit
    if tstamp > timelimit
      s.termination = "Time Limit"
      return s
    end
  end

  s.termination = "Max Iterations"
  return s
end

function forwardstep(graph::PlasmoGraph, cuts::Array{Symbol,1}, updatebound::Bool)
  levels = graph.attributes[:levels]
  numlevels = length(levels)
  for level in 1:numlevels
    nodeslevel = levels[level]
    for node in nodeslevel
      nodeworkflow(node,graph,cuts,updatebound)
    end
  end
  LB = graph.attributes[:LB]
  if updatebound
    UB = sum(node.attributes[:preobjval] for node in values(graph.nodes))
    graph.attributes[:UB] = UB
  else
    UB = graph.attributes[:UB]
  end
  return LB,UB
end

function nodeworkflow(node::PlasmoNode, graph::PlasmoGraph, cuts::Array{Symbol,1}, updatebound::Bool)
  # 1. Add cuts
  generatecuts(node,graph)
  # 2. Take x
  takex(node)
  # 3. solve
  if :LP in cuts
    solvelprelaxation(node)
  end
  if :Root in cuts
    solverootrelaxation(node)
  end
  if updatebound
    solvenodemodel(node,graph)
  end
  # 4. put x
  putx(node,graph)
  # 5. put cuts and nodebound
  putcutdata(node,graph,cuts)
end

function solvelprelaxation(node::PlasmoNode)
  model = getmodel(node)
  status = solve(model, relaxation = true)

  @assert status == :Optimal

  dualconstraints = node.attributes[:linkconstraints]

  λnode = getdual(dualconstraints)
  nodebound = getobjectivevalue(model)

  node.attributes[:bound] = nodebound
  node.attributes[:λ] = λnode

  return status
end

function solverootrelaxation(node::PlasmoNode)
  sp = getmodel(node)
  if length(sp.linconstrDuals) == 0
    solve(sp, relaxation=true)
  end
  lpfile = joinpath(tmpdir,"nodemodel.lp")
  writeLP(sp,lpfile)
  run(`cpxgetroot $lpfile 0 1`)
  lp = Model(solver=CPLEX.CplexSolver(CPX_PARAM_PREIND=0))
  lp.internalModel = MathProgBase.LinearQuadraticModel(lp.solver)
  MathProgBase.loadproblem!(lp.internalModel,"node0.lp")
  MathProgBase.optimize!(lp.internalModel)

  run(`mv node0.lp $tmpdir/`)

  dualconstraints = node.attributes[:linkconstraints]

  rootduals = MathProgBase.getconstrduals(lp.internalModel)
  sp.linconstrDuals = MathProgBase.getconstrduals(lp.internalModel)[1:length(sp.linconstrDuals)]

  λnode = getdual(dualconstraints)
  nodebound = MathProgBase.getobjval(lp.internalModel)

  node.attributes[:bound] = nodebound
  node.attributes[:λ] = λnode
end

function solvenodemodel(node::PlasmoNode,graph::PlasmoGraph)
  model = getmodel(node)
  solve(model)
  if in_degree(graph,node) == 0 # Root node
    graph.attributes[:LB] = getobjectivevalue(model)
  end
  node.attributes[:preobjval] = getvalue(model.ext[:preobj])
end

function takex(node::PlasmoNode)
  xinvals = node.attributes[:xin]
  xinvars = node.attributes[:xinvars]
  if length(xinvals) > 0
    fix.(xinvars,xinvals)
  end
end

function putx(node::PlasmoNode,graph::PlasmoGraph)
  childvars = node.attributes[:childvars]
  children = out_neighbors(graph,node)
  length(children) == 0 && return true

  for child in children
    xnode = getvalue(childvars[getnodeindex(graph,child)])
    child.attributes[:xin] = xnode
  end
end

function putcutdata(node::PlasmoNode,graph::PlasmoGraph,cuts::Array{Symbol,1})
  parents = in_neighbors(graph,node)
  length(parents) == 0 && return true
  parent = parents[1]    # Assume only one parent
  parentcuts = parent.attributes[:cutdata]
  θk = node.attributes[:bound]
  λk = node.attributes[:λ]
  xk = node.attributes[:xin]
  nodeindex = getnodeindex(graph,node)
  if :LP in cuts || :Root in cuts
    bcd = BendersCutData(θk, λk, xk)
    push!(parentcuts[nodeindex],bcd)
  end
  if :LLinteger in cuts
    llintcd = LLIntegerCutData(θlb,xk)
    push!(parentcuts[nodeindex],llintcd)
  end
  if :Integer in cuts
    intcd = IntegerCutData(xk)
    push!(parentcuts[nodeindex],intcd)
  end
end

function generatecuts(node::PlasmoNode,graph::PlasmoGraph)
  children = out_neighbors(graph,node)
  length(children) == 0 && return true

  cutdataarray = node.attributes[:cutdata]
  for child in children
    childindex = getnodeindex(graph,child)
    while length(cutdataarray[childindex]) > 0
      cutdata = pop!(cutdataarray[childindex])
      if typeof(cutdata) == BendersCutData
        generatebenderscut(node,cutdata,childindex)
      elseif typeof(cutdata) == LLIntegerCutData
        generateLLintegercut(node,cutdata)
      elseif typeof(cutdata) == IntegerCutData
        generateintegercut(node,cutdata)
      end
    end
  end
end

function generatebenderscut(node::PlasmoNode, cd::BendersCutData,index)
  model = getmodel(node)
  θ = getindex(model, :θ)
  x = node.attributes[:childvars][index]
  @constraint(model, θ[index] >= cd.θk + cd.λk'*(x - cd.xk))
end


function identifylevels(graph::Plasmo.PlasmoGraph)
  #Create lists of root and leaf nodes in graph
  roots = graph.attributes[:roots] = []
  leaves = graph.attributes[:leaves] = []
  #Create dictionary to keep track of levels of nodes
  levels = graph.attributes[:levels] = Dict()
  #Iterate through every node to check for root/leaf nodes
  for node in values(graph.nodes)
    node.attributes[:xin] = []
    node.attributes[:λ] = []
    node.attributes[:bound] = NaN
    node.attributes[:xinvars] = []
    node.attributes[:preobjval] = NaN
    node.attributes[:linkconstraints] = []
    #If the node does not have parents it is a root node
    if in_degree(graph,node) == 0
      push!(roots,node)
    end
    #If the node does not have children it is a leaf node
    if out_degree(graph,node) == 0
      push!(leaves,node)
    end
  end

  #Start mapping level from the root nodes
  current = roots
  level = 1
  levels[1]  = roots
  while length(current) > 0
    levels[level] = current
    children = []
    for node in current
        push!(children,out_neighbors(graph,node)...)
        node.attributes[:childvars] = Dict(getnodeindex(graph,child) => [] for child in out_neighbors(graph,node))
        node.attributes[:cutdata] = Dict(getnodeindex(graph,child) => CutData[] for child in out_neighbors(graph,node))
    end
    current = children
    level += 1
  end
  graph.attributes[:numlevels] = level - 1
end

function bdprepare(graph::Plasmo.PlasmoGraph)
  if haskey(graph.attributes,:preprocessed)
    return true
  end

  identifylevels(graph)
  graph.attributes[:normalized] = normalizegraph(graph)
  graph.attributes[:mflat] = create_flat_graph_model(graph)
  setsolver(graph.attributes[:mflat],graph.solver)

  links = getlinkconstraints(graph)
  numlinks = length(links)

  for index in 1:length(graph.nodes)
    node = graph.nodes[index]
    model = getmodel(node)
    if model.solver == JuMP.UnsetSolver()
      model.solver = graph.solver
    end
    model.ext[:preobj] = model.obj
    #Add theta to parent nodes
    if out_degree(graph,node) != 0
      childrenindices = [getnodeindex(graph,child) for child in out_neighbors(graph,node)]
      sort!(childrenindices)
      @variable(model, θ[i in childrenindices] >= -1e6)
      model.obj += sum(θ[i] for i in childrenindices)
    end
  end
  #Add dual constraint to child nodes using the linking constraints
  for (numlink,link) in enumerate(links)
    #Take the two variables of the constraint
    var1 = link.terms.vars[1]
    var2 = link.terms.vars[2]
    #Determine which nodes they belong to
    nodeV1 = getnode(var1)
    nodeV2 = getnode(var2)
    #Set the order of the nodes
    if ischildnode(graph,nodeV1,nodeV2)
      childnode = nodeV1
      childvar = var1
      parentnode = nodeV2
      parentvar = var2
    else
      childnode = nodeV2
      childvar = var2
      parentnode = nodeV1
      parentvar = var1
    end
    childindex = getnodeindex(graph,childnode)
    childmodel = getmodel(childnode)
    push!(parentnode.attributes[:childvars][childindex],parentvar)
    valbar = @variable(childmodel)
    setname(valbar,"varlbar$numlink")
    push!(childnode.attributes[:xinvars],valbar)
    conref = @constraint(childmodel, valbar - childvar == 0)
    push!(childnode.attributes[:linkconstraints], conref)
  end
  graph.attributes[:preprocessed] = true
end

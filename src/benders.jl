
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

(==)(cd1::BendersCutData,cd2::BendersCutData) = (cd1.θk == cd2.θk) && (cd1.λk == cd2.λk) && (cd1.xk == cd2.xk)
(==)(cd1::LLIntegerCutData,cd2::LLIntegerCutData) = (cd1.θlb == cd2.θlb) &&  (cd1.yk == cd2.yk)
(==)(cd1::IntegerCutData,cd2::IntegerCutData) = (cd1.yk == cd2.yk)

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
  update!(graph.attributes[:LB],LB)
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

    tstamp = time() - starttime
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

    if fetch(graph.attributes[:stalled])
      s.termination = "Stalled"
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
      solveprimalnode(node,cuts,updatebound)
    end
  end
  LB = fetch(graph.attributes[:LB])
  if updatebound
    UB = sum(node.attributes[:preobjval] for node in values(graph.nodes))
    graph.attributes[:UB] = UB
  else
    UB = graph.attributes[:UB]
  end
  return LB,UB
end

function solveprimalnode(node::PlasmoNode, cuts::Array{Symbol,1}, updatebound::Bool)
  # 1. Add cuts
  generatecuts(node)
  # 2. Take x
  takex(node)
  # 3. solve
  if :LP in cuts
    solvelprelaxation(node)
  end
  if updatebound
    solvenodemodel(node)
  end
  # 4. put x
  putx(node)
  # 5. put cuts and nodebound
  putcutdata(node,cuts)
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


function solvenodemodel(node::PlasmoNode)
  model = getmodel(node)
  solve(model)
  if node.attributes[:numparents] == 0 # Root node
     update!(node.attributes[:graphLB], getobjectivevalue(model))
  end
  node.attributes[:preobjval] = getvalue(model.ext[:preobj])
end

function takex(node::PlasmoNode)
  node.attributes[:numparents] == 0 && return true
  xinvals = fetch(node.attributes[:xin])
  xinvars = node.attributes[:xinvars]
  if length(xinvals) > 0
    fix.(xinvars,xinvals)
  end
end

function putx(node::PlasmoNode)
  node.attributes[:numchildren] == 0 && return true
  childvars = node.attributes[:childvars]
  for childindex in keys(childvars)
    xnode = getvalue(childvars[childindex])
    put!(node.attributes[:xout][childindex],xnode)
  end
end

function putcutdata(node::PlasmoNode,cuts::Array{Symbol,1})
  node.attributes[:numparents] == 0 && return true
  θk = node.attributes[:bound]
  λk = node.attributes[:λ]
  xk = take!(node.attributes[:xin])
  if :LP in cuts || :Root in cuts
    bcd = BendersCutData(θk, λk, xk)
    put!(node.attributes[:cutdataout],bcd)
  end
  if :LLinteger in cuts
    llintcd = LLIntegerCutData(θlb,xk)
    put!(node.attributes[:cutdataout],llintcd)
  end
  if :Integer in cuts
    intcd = IntegerCutData(xk)
    put!(node.attributes[:cutdataout],intcd)
  end
end

function generatecuts(node::PlasmoNode)
  node.attributes[:numchildren] == 0 && return true

  cutdataarray = node.attributes[:cutdata]
  previouscuts = node.attributes[:prevcuts]
  thisitercuts = Dict()
  samecuts = Dict()
  for childindex in keys(cutdataarray)
    thisitercuts[childindex] = CutData[]
    samecuts[childindex] = Bool[]
    while isready(cutdataarray[childindex]) 
      cutdata = take!(cutdataarray[childindex])
      samecut = in(cutdata,previouscuts[childindex])
      push!(samecuts[childindex],samecut)
      samecut && continue
      if typeof(cutdata) == BendersCutData
        generatebenderscut(node,cutdata,childindex)
      elseif typeof(cutdata) == LLIntegerCutData
        generateLLintegercut(node,cutdata)
      elseif typeof(cutdata) == IntegerCutData
        generateintegercut(node,cutdata)
      end
      push!(thisitercuts[childindex],cutdata)
    end
    samecuts[childindex] = reduce(*,samecuts[childindex]) && length(samecuts[childindex]) > 0
  end
  node.attributes[:prevcuts] = thisitercuts
  nodesamecuts = collect(values(samecuts))
    stalled = reduce(*,nodesamecuts)
    node.attributes[:stalled] = stalled
  stalled && warn("Root node stalled")
  if node.attributes[:numparents] == 0 && stalled     
    update!(node.attributes[:graphstalled],true)
  end
end

function generatebenderscut(node::PlasmoNode, cd::BendersCutData,index)
  model = getmodel(node)
  θ = getindex(model, :θ)
  x = node.attributes[:childvars][index]
  @constraint(model, θ[index] >= cd.θk + cd.λk'*(cd.xk - x))
end


function identifylevels(graph::Plasmo.PlasmoGraph)
  #Create lists of root and leaf nodes in graph
  roots = graph.attributes[:roots] = []
  leaves = graph.attributes[:leaves] = []
  #Create dictionary to keep track of levels of nodes
  levels = graph.attributes[:levels] = Dict()
  #Iterate through every node to check for root/leaf nodes
  for node in values(graph.nodes)
    node.attributes[:λ] = []
    node.attributes[:bound] = NaN
    node.attributes[:xinvars] = []
    node.attributes[:preobjval] = NaN
    node.attributes[:linkconstraints] = []
    #If the node does not have parents it is a root node
      if in_degree(graph,node) == 0
         LBchannel = RemoteChannel(()->Channel{Float64}(1))
         put!(LBchannel,-Inf)
         graph.attributes[:LB] = LBchannel
         node.attributes[:numparents] = 0
         node.attributes[:graphLB] = LBchannel
         stalledchannel = RemoteChannel(()->Channel{Bool}(1))
         put!(stalledchannel,false)
         node.attributes[:graphstalled] = stalledchannel
         graph.attributes[:stalled] = stalledchannel
         push!(roots,node)
      end
    #If the node does not have children it is a leaf node
      if out_degree(graph,node) == 0
          node.attributes[:numchildren] = 0
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
        node.attributes[:childvars] = Dict()
        node.attributes[:xout] = Dict()
        node.attributes[:cutdata] = Dict()
        node.attributes[:prevcuts] = Dict()
        node.attributes[:numchildren] = 0
        for child in children
            childindex = getindex(graph,child)
            node.attributes[:childvars][childindex] = []

            xchannel = RemoteChannel(()->Channel{Array{Float64,1}}(1))
            node.attributes[:xout][childindex] = xchannel
            child.attributes[:xin] = xchannel

            cutchannel = RemoteChannel(()->Channel{CutData}(10))
            node.attributes[:cutdata][childindex] = cutchannel
            child.attributes[:cutdataout] = cutchannel

            node.attributes[:childvars][childindex] = []
            node.attributes[:numchildren] += 1
            child.attributes[:numparents] = 1
        end
        node.attributes[:stalled] = false
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

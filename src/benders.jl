include("dualliftproject.jl")
"""
bendersolve
"""
function bendersolve(graph::ModelGraph; max_iterations::Int64=10, cuts::Array{Symbol,1}=[:LP], ϵ=1e-5,UBupdatefrequency=1,timelimit=3600,verbose=false)
  starttime = time()
  s = Solution(method=:benders)
  setattribute(graph, :solution, s)
  updatebound = true

  verbose && info("Preparing graph")
  bdprepare(graph, cuts)
  
  n = getattribute(graph, :normalized)

  verbose && info("Solve relaxation and set LB")
  mf = getattribute(graph, :mflat)
  solve(mf,relaxation=true)
  LB = getobjectivevalue(getattribute(graph, :mflat))
  UB = Inf

  # Set bound to root node
  rootnode = getattribute(graph, :roots)[1]
  rootmodel = getmodel(rootnode)
  @constraint(rootmodel, rootmodel.obj.aff >= LB)

  # Begin iterations
  for i in 1:max_iterations
    tic()
    updatebound = ((i-1) % UBupdatefrequency) == 0
    LB,UB = forwardstep(graph, cuts, updatebound)

    tstamp = time() - starttime

    itertime = toc()
    if n == 1
      saveiteration(s,tstamp,[UB,LB,itertime,tstamp],n)
    else
      saveiteration(s,tstamp,[n*LB,n*UB,itertime,tstamp],n)
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

    if getattribute(graph, :stalled)
      s.termination = "Stalled"
      return s
    end
  end

  s.termination = "Max Iterations"
  return s
end

function forwardstep(graph::ModelGraph, cuts::Array{Symbol,1}, updatebound::Bool)
  levels = getattribute(graph, :levels)
  numlevels = length(levels)
  for level in 1:numlevels
    nodeslevel = levels[level]
    for node in nodeslevel
      solveprimalnode(node,graph,cuts,updatebound)
    end
  end
  LB = getattribute(graph, :LB)
  if updatebound
    iterUB = sum(getattribute(node, :preobjval) for node in getnodes(graph))
    setattribute(graph, :iterUB, iterUB)
    UB = min(getattribute(graph, :UB),iterUB)
    setattribute(graph, :UB, UB)
  else
    UB = getattribute(graph, :UB)
  end
  return LB,UB
end

function solveprimalnode(node::ModelNode, graph::ModelGraph, cuts::Array{Symbol,1}, updatebound::Bool)
  # 1. Add cuts
  generatecuts(node,graph)
  # 2. Take x
  takex(node)
  # 3. solve
  if :LP in cuts
    solvelprelaxation(node)
  end
  if :LIFT in cuts && in_degree(graph, node) != 0
    solveliftandprojectrelaxation(node, graph)
  end
  if updatebound
    solvenodemodel(node,graph)
  end
  # 4. put x
  putx(node,graph)
  # 5. put cuts and nodebound
  putcutdata(node,graph,cuts)
end

function solvelprelaxation(node::ModelNode)
  model = getmodel(node)
  status = solve(model, relaxation = true)

  @assert status == :Optimal

  dualconstraints = getattribute(node, :linkconstraints)

  λnode = getdual(dualconstraints)
  nodebound = getobjectivevalue(model)

  setattribute(node, :bound, nodebound)
  setattribute(node, :λ, λnode)

  return status
end

function solvenodemodel(node::ModelNode,graph::ModelGraph)
  model = getmodel(node)
  solve(model)
  if in_degree(graph,node) == 0 # Root node
    setattribute(graph, :LB, getobjectivevalue(model))
  end
  setattribute(node, :preobjval, JuMP.getvalue(model.ext[:preobj]))
end

function takex(node::ModelNode)
  xinvals = getattribute(node, :xin)
  xinvars = getattribute(node, :xinvars)
  if length(xinvals) > 0
    fix.(xinvars,xinvals)
  end
end

function putx(node::ModelNode,graph::ModelGraph)
  childvars = getattribute(node,:childvars)
  children = out_neighbors(graph,node)
  length(children) == 0 && return true

  for child in children
    xnode = JuMP.getvalue(childvars[getindex(graph,child)])
    #round xnode to bounds (sometimes numerical errors can occur)
    # for i in 1:length(childvars[getnodeindex(graph,child)])
    #   var = childvars[getnodeindex(graph,child)][i]
    #   ub = getupperbound(var)
    #   lb = getlowerbound(var)
    #   category = getcategory(var)
    #   if xnode[i] > ub
    #     xnode[i] = ub
    #   end
    #   if xnode[i] < lb
    #     xnode[i] = lb
    #   end
    #   if category == :Bin
    #     if xnode[i] < 1e-4
    #       xnode[i] = 0
    #     else
    #       xnode[i] = 1
    #     end
    #   end
    # end
    setattribute(child,:xin, xnode)
  end
end

function putcutdata(node::ModelNode,graph::ModelGraph,cuts::Array{Symbol,1})
  parents = in_neighbors(graph,node)
  length(parents) == 0 && return true
  parent = parents[1]    # Assume only one parent
  parentcuts = getattribute(parent, :cutdata)
  θk = getattribute(node,:bound)
  λk = getattribute(node,:λ)
  xk = getattribute(node,:xin)
  nodeindex = getindex(graph,node)
  if :LP in cuts || :Root in cuts || :GMI in cuts || :LIFT in cuts
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

function generatecuts(node::ModelNode,graph::ModelGraph)
  children = out_neighbors(graph,node)
  length(children) == 0 && return true

  cutdataarray = getattribute(node,:cutdata)
  previouscuts = getattribute(node,:prevcuts)
  thisitercuts = Dict()
  samecuts = Dict()
  for child in children
    childindex = getindex(graph,child)
    thisitercuts[childindex] = CutData[]
    samecuts[childindex] = Bool[]

    while length(cutdataarray[childindex]) > 0
      cutdata = pop!(cutdataarray[childindex])
      samecut = in(cutdata,previouscuts[childindex])
      push!(samecuts[childindex],samecut)
      samecut && continue
      if typeof(cutdata) == BendersCutData
          generatebenderscut(node, cutdata, childindex)
      elseif typeof(cutdata) == LLIntegerCutData
          generateLLintegercut(node,cutdata)
      elseif typeof(cutdata) == IntegerCutData
          generateintegercut(node,cutdata)
      elseif typeof(cutdata) == LagrangeCrossCutData
          generatelagrangecrosscut(node, cutdata, childindex)
      end
      push!(thisitercuts[childindex],cutdata)
    end
    samecuts[childindex] = reduce(*,samecuts[childindex]) && length(samecuts[childindex]) > 0
  end
  setattribute(node,:prevcuts, thisitercuts)
  nodesamecuts = collect(values(samecuts))
  cuts = getattribute(graph, :cuts)
  if :LP in cuts
    setattribute(node,:stalled, reduce(*,nodesamecuts))
    #bound does not improve for 5 iterations is also considered as stalled
    if length(getattribute(graph, :solution).iterbound) > 6 && abs(getattribute(graph, :LB) - getattribute(graph, :solution).iterbound[end-5]) < 1e-3 && abs(getattribute(graph, :solution).iterval[end] - getattribute(graph, :solution).iterval[end-5]) < 1e-3
      setattribute(node, :stalled, true)
    end
    getattribute(node, :stalled) && warn("Node stalled")
  end
#only stalled when LP has already stalled
  if (:GMI in cuts || :LIFT in cuts) && getattribute(node, :LP_stalled)
    setattribute(node, :stalled, reduce(*,nodesamecuts))
    if length(getattribute(graph, :solution).iterbound) > 6 && abs(getattribute(graph, :LB) - getattribute(graph, :solution).iterbound[end-5]) < 1e-3 && abs(getattribute(graph, :solution).iterval[end] - getattribute(graph, :solution).iterval[end-5]) < 1e-3 && getattribute(node, :LP_stalled_iterations) > 2
      setattribute(node, :stalled, true)
    end
    setattribute(node, :LP_stalled_iterations, getattribute(node, :LP_stalled_iterations) + 1 )
    getattribute(node, :stalled) && warn("Node  stalled")
  end
  if (:GMI in cuts || :LIFT in cuts) && (!getattribute(node, :LP_stalled))
    setattribute(node, :LP_stalled, reduce(*,nodesamecuts))
    if length(getattribute(graph, :solution).iterbound) > 6 && abs(getattribute(graph, :LB) - getattribute(graph, :solution).iterbound[end-5]) < 1e-3 && abs(getattribute(graph, :solution).iterval[end] - getattribute(graph, :solution).iterval[end-5]) < 1e-3
      setattribute(node, :LP_stalled, true)
    end
    println("LP_stalled status")
    println(getattribute(node, :LP_stalled))
    if getattribute(node, :LP_stalled)
      setattribute(node, :LP_stalled_iterations, 1)
    end
  end
  if in(node,getattribute(graph, :roots)) && getattribute(node, :stalled)
    setattribute(graph, :stalled, true)
  end


end

function generatebenderscut(node::ModelNode, cd::BendersCutData,index)
  model = getmodel(node)
  θ = getindex(model, :θ)
  x = getattribute(node, :childvars)[index]
  @constraint(model, θ[index] >= cd.θk + cd.λk'*(cd.xk - x))
end

function generatelagrangecrosscut(node, cd::LagrangeCrossCutData, index)
    model = getmodel(node)
    θ = getindex(model, :θ)
    x = getattribute(node, :childvars)[index]
    @constraint(model, θ[index] >= cd.zk + cd.λk'*x)
end

function identifylevels(graph::ModelGraph)
  #Create lists of root and leaf nodes in graph
  setattribute(graph, :roots, ModelNode[])
  roots = getattribute(graph, :roots)
  setattribute(graph, :leaves, ModelNode[])
  leaves = getattribute(graph, :leaves)
  #Create dictionary to keep track of levels of nodes
  setattribute(graph, :levels, Dict())
  levels = getattribute(graph, :levels)
  #Iterate through every node to check for root/leaf nodes
  for node in getnodes(graph)
    setattribute(node,:xin, [])
    setattribute(node,:λ, [])
    setattribute(node,:bound, NaN)
    setattribute(node,:xinvars, [])
    setattribute(node,:preobjval, NaN)
    setattribute(node,:linkconstraints, [])
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
        setattribute(node,:childvars, Dict(getindex(graph,child) => [] for child in out_neighbors(graph,node)))
        setattribute(node,:cutdata, Dict(getindex(graph,child) => CutData[] for child in out_neighbors(graph,node)))
        setattribute(node,:prevcuts, Dict(getindex(graph,child) => CutData[] for child in out_neighbors(graph,node)))
        setattribute(node,:stalled, false)
    end
    current = children
    level += 1
  end
  setattribute(graph,:numlevels,level - 1)
end

function bdprepare(graph::ModelGraph, cuts::Array{Symbol,1}=[:LP])
  if Plasmo.hasattribute(graph,:preprocessed)
    return true
  end

  identifylevels(graph)
  setattribute(graph, :normalized, normalizegraph(graph))
  setattribute(graph, :stalled, false)
  setattribute(graph, :mflat ,create_jump_graph_model(graph))
  setattribute(graph, :UB, Inf)
  setattribute(graph, :cuts, cuts)
  JuMP.setsolver(getattribute(graph, :mflat) ,getsolver(graph))

  links = getlinkconstraints(graph)
  numlinks = length(links)

  for index in 1:length(getnodes(graph))
    node = getnode(graph, index)
    model = getmodel(node)
    if model.solver == JuMP.UnsetSolver()
      model.solver = getsolver(graph)
    end

#add code for :LIFT========

    #if :GMI in cuts || :LIFT in cuts
    #create upper bound for variables if GMI or LIFT is in cuts
      for col in 1:length(model.colUpper)
        if model.colUpper[col] > 1e10
          setupperbound(Variable(model, col), 1e10)

        end
      end
      #add attributes to check if LP cuts has stalled
    #end
    setattribute(node, :LP_stalled, false)

    #create standard matrix for lift and project cuts
    if :LIFT in cuts && in_degree(graph, node) != 0
      getmatrixform(node)
    end
#finish add code =========

    model.ext[:preobj] = model.obj
    setattribute(node, :cgmodel, deepcopy(model))
    #Add theta to parent nodes
    if out_degree(graph,node) != 0
      childrenindices = [getindex(graph,child) for child in out_neighbors(graph,node)]
      sort!(childrenindices)
      @variable(model, θ[i in childrenindices] >= -1e6)
      model.obj += sum(θ[i] for i in childrenindices)
    end
    #get the index of linking variables
    if in_degree(graph, node) != 0
      setattribute(node, :linking_vars_indices, [])
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
    childindex = getindex(graph,childnode)
    childmodel = getmodel(childnode)
    push!(getattribute(parentnode, :childvars)[childindex],parentvar)
    linkvar = @variable(childmodel)
    setname(linkvar,"linkvar$numlink")
    push!(getattribute(childnode, :xinvars),linkvar)
    conref = @constraint(childmodel, linkvar - childvar == 0)
    push!(getattribute(childnode, :linkconstraints), conref)
    #store child var index
    push!(getattribute(childnode, :linking_vars_indices), childvar.col)
    push!(getattribute(childnode, :linking_vars_indices), linkvar.col)
  end

  #set up CGLPs
  if :LIFT in cuts
    for node in values(getnodes(graph))
      if in_degree(graph, node) != 0
        setCGLP(node, graph)
      end
    end
  end

  setattribute(graph, :preprocessed, true)
end

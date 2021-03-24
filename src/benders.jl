"""
bendersolve
"""
function bendersoptimize!(graph::OptiGraph; 
  max_iterations::Int64=10, 
  cuts::Array{Symbol,1}=[:LP], 
  ϵ=1e-5,UBupdatefrequency=1,
  timelimit=3600, 
  singleline=true,
  verbose=true)

  starttime = time()
  global tmpdir = "/tmp/RootNode" # mktempdir()
  s = Solution(method=:benders)
  updatebound = true

  bdprepare(graph)
  n = getattribute(graph,:normalized)
  LB = initialrelaxation(graph)
  update!(getattribute(graph,:LB),LB)
  
  # Set bound to root node
  rootnode = getattribute(graph, :roots)[1]
  rootmodel = getmodel(rootnode)
  @constraint(rootmodel, objective_function(rootmodel) >= LB)

  # Send nodes to workers
  nprocs() > 1 && sendnodestoworkers(graph)
  # Begin iterations
  for i in 1:max_iterations
    tic = time()
    updatebound = ((i-1) % UBupdatefrequency) == 0
    LB,UB = forwardstep(graph, cuts, updatebound)

    tstamp = time() - starttime
    if n == 1
      saveiteration(s,tstamp,[UB,LB,time() - tic,tstamp],n)
    else
      saveiteration(s,tstamp,[n*LB,n*UB,time() - tic,tstamp],n)
    end
    if verbose 
      i == 1 && singleline && printheader(graph, :Benders)
      printiterationsummary(s,singleline=singleline)
    end

    if abs(UB-LB) < ϵ
      s.termination = "Optimal"
      return s
    end

    # Check time limit
    if tstamp > timelimit
      s.termination = "Time Limit"
      return s
    end

    if fetch(graph.obj_dict[:stalled])
      s.termination = "Stalled"
      return s
    end
  end

  s.termination = "Max Iterations"
  return s
end

function forwardstep(graph::OptiGraph, cuts::Array{Symbol,1}, updatebound::Bool)
  levels = graph.obj_dict[:levels]
  numlevels = length(levels)
  for level in 1:numlevels
    nodeslevel = levels[level]
    worker = repeat(workers(),outer=Int(ceil(length(nodeslevel)/nworkers())))
    @sync for (k,node) in enumerate(nodeslevel)
      remotecall_fetch(solveprimalnode,worker[k],node,cuts,updatebound)
    end
  end
  LB = fetch(graph.obj_dict[:LB])
  if updatebound
    UB = sum(fetch(graph.obj_dict[:preobjval][nodelabel]) for nodelabel in keys(graph.obj_dict[:preobjval]))
    graph.obj_dict[:UB] = UB
  else
    UB = graph.obj_dict[:UB]
  end
  return LB,UB
end

function solveprimalnode(node::OptiNode, cuts::Array{Symbol,1}, updatebound::Bool)
  # 1. Add cuts
  generatecuts(node)
  # 2. Take x
  takex(node)
  if node[:updated]
    # 3. solve
    if :LP in cuts
      solvelprelaxation(node)
    end
    if node[:primalstatus] == MOI.FEASIBLE_POINT
      if updatebound
        solvenodemodel(node)
      end
      # 4. put x
      putx(node)
    end
    # 5. put cuts and nodebound
    putcutdata(node,cuts)
  else
    # Flush x_in
    xin = take!(node[:xin])
  end
end


function solvelprelaxation(node::OptiNode)
  model = getmodel(node)
  unrelax = relax_integrality(model)
  optimize!(model)

  node[:primalstatus] = primal_status(model)
  node[:dualstatus] = dual_status(model)

  node[:primalstatus] == MOI.NO_SOLUTION && node[:dualstatus] == MOI.NO_SOLUTION && error("Infeasible model at node, solver does not provide infeasibility certificate")
  dualconstraints = node[:linkconstraints]

  λnode = dual.(dualconstraints)
  nodebound = node[:primalstatus] == MOI.FEASIBLE_POINT ? objective_value(model) : 0

  node[:bound] = nodebound
  node[:λ] = λnode

  unrelax()
end


function solvenodemodel(node::OptiNode)
  model = getmodel(node)
  optimize!(model)
  if node[:numparents] == 0 # Root node
     update!(node[:graphLB], objective_value(model))
  end
  update!(node[:preobjval],getvalue(model.ext[:preobj]))
end


function takex(node::OptiNode)
  node[:numparents] == 0 && return true
  xinvals = fetch(node[:xin])
  if xinvals == node[:prevx] && node[:cutsgenerated] == 0
    node[:updated] = false
  else
    xinvars = node[:xinvars]
    if length(xinvals) > 0
      fix.(xinvars,xinvals)
    end
    node[:updated] = true
    node[:prevx] = xinvals
  end
end

function putx(node::OptiNode)
  node[:numchildren] == 0 && return true
  childvars = node[:childvars]
  for childindex in keys(childvars)
    xnode = value.(childvars[childindex])
    put!(node[:xout][childindex],xnode)
  end
end

function putcutdata(node::OptiNode,cuts::Array{Symbol,1})
  node[:numparents] == 0 && return true
  θk = node[:bound]
  λk = node[:λ]
  xk = take!(node[:xin])
  if node[:dualstatus] == MOI.INFEASIBILITY_CERTIFICATE
    fcd = FeasibilityCutData(λk, xk)
    put!(node[:cutdataout],fcd)
  else
    if :LP in cuts || :Root in cuts
      bcd = BendersCutData(θk, λk, xk)
      put!(node[:cutdataout],bcd)
    end
    if :LLinteger in cuts
      llintcd = LLIntegerCutData(θlb,xk)
      put!(node[:cutdataout],llintcd)
    end
    if :Integer in cuts
      intcd = IntegerCutData(xk)
      put!(node[:cutdataout],intcd)
    end
  end
end

function generatecuts(node::OptiNode)
  node[:numchildren] == 0 && return true
  node[:cutsgenerated] = 0
  cutdataarray = node[:cutdata]
  previouscuts = node[:prevcuts]
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
      elseif typeof(cutdata) == FeasibilityCutData
        generatefeasibilitycut(node,cutdata,childindex)
      elseif typeof(cutdata) == LLIntegerCutData
        generateLLintegercut(node,cutdata)
      elseif typeof(cutdata) == IntegerCutData
        generateintegercut(node,cutdata)
      elseif typeof(cutdata) == LagrangeCrossCutData
        generatelagrangecrosscut(node, cutdata, childindex)
      end
      node[:cutsgenerated] += 1
      push!(thisitercuts[childindex],cutdata)
    end
    samecuts[childindex] = reduce(*,samecuts[childindex]) && length(samecuts[childindex]) > 0
  end
  node[:prevcuts] = thisitercuts
  nodesamecuts = collect(values(samecuts))
    stalled = reduce(*,nodesamecuts)
    node[:stalled] = stalled
  stalled && warn("Root node stalled")
  if node[:numparents] == 0 && stalled     
    update!(node[:graphstalled],true)
  end
end


function generatebenderscut(node::OptiNode, cd::BendersCutData,index)
  model = getmodel(node)
  θ = model[:θ]
  x = node[:childvars][index]
  @constraint(model, θ[index] >= cd.θk + sum(cd.λk[i]*(cd.xk[i] - x[i]) for i in keys(x)))
end


function generatefeasibilitycut(node::OptiNode, cd::FeasibilityCutData,index)
  model = getmodel(node)
  θ = model[:θ]
  x = node[:childvars][index]
  @constraint(model, θ[index] >= sum(cd.λk[i]*(cd.xk[i] - x[i]) for i in keys(x)))
end


function generatelagrangecrosscut(node::OptiNode, cd::LagrangeCrossCutData, index)
  model = getmodel(node)
  θ = model[:θ]
  x = getattribute(node, :childvars)[index]
  @constraint(model, θ[index] >= cd.zk + sum(cd.λk[i]*x[i] for i in keys(x)))
end


function identifylevels(graph::OptiGraph)
  #Create lists of root and leaf nodes in graph
  roots = graph.obj_dict[:roots] = []
  leaves = graph.obj_dict[:leaves] = []
  #Create dictionary to keep track of levels of nodes
    levels = graph.obj_dict[:levels] = Dict()
  #Create channel to receive preobj
  graph.obj_dict[:preobjval] = Dict(node.label => RemoteChannel(()->Channel{Float64}(1)) for node in values(getnodes(graph)))
  #Iterate through every node to check for root/leaf nodes
  for node in getnodes(graph)
    node[:λ] = []
    node[:bound] = NaN
    node[:xinvars] = []
    node[:preobjval] = graph.obj_dict[:preobjval][node.label]
    put!(node[:preobjval],NaN)  
    node[:linkconstraints] = []
    #If the node does not have parents it is a root node
    if in_degree(graph ,node) == 0
        LBchannel = RemoteChannel(()->Channel{Float64}(1))
        put!(LBchannel,-Inf)
        graph.obj_dict[:LB] = LBchannel
        node[:numparents] = 0
        node[:graphLB] = LBchannel
        stalledchannel = RemoteChannel(()->Channel{Bool}(1))
        put!(stalledchannel,false)
        node[:graphstalled] = stalledchannel
        graph.obj_dict[:stalled] = stalledchannel
        push!(roots, node)
    end
    #If the node does not have children it is a leaf node
    if out_degree(graph, node) == 0
        node[:numchildren] = 0
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
        push!(children, out_neighbors(graph,node)...)
        node[:childvars] = Dict()
        node[:xout] = Dict()
        node[:cutdata] = Dict()
        node[:prevcuts] = Dict()
        node[:numchildren] = 0
        node[:stalled] = false
        node[:cutsgenerated] = 0
        node[:updated] = true
        node[:prevx] = []
        node[:primalstatus] = MOI.FEASIBLE_POINT
        node[:dualstatus] = MOI.FEASIBLE_POINT
        
        for child in out_neighbors(graph,node)
            childindex = graph[child]
            node[:childvars][childindex] = []

            xchannel = RemoteChannel(()->Channel{Array{Float64,1}}(1))
            node[:xout][childindex] = xchannel
            child[:xin] = xchannel

            cutchannel = RemoteChannel(()->Channel{CutData}(10))
            node[:cutdata][childindex] = cutchannel
            node[:prevcuts][childindex] = []
            child[:cutdataout] = cutchannel

            node[:childvars][childindex] = []
            node[:numchildren] += 1
            child[:numparents] = 1
        end
    end
    current = children
    level += 1
  end
  graph.obj_dict[:numlevels] = level - 1
end

function bdprepare(graph::OptiGraph)
  if haskey(graph.obj_dict,:preprocessed)
    return true
  end

  identifylevels(graph)
  graph.obj_dict[:normalized] = normalizegraph(graph)
  optinode,  = combine(graph)
  graph.obj_dict[:mflat] = getmodel(optinode)
  graph.obj_dict[:numlinks] = length(getlinkconstraints(graph))
  graph.obj_dict[:λ] = [zeros(graph.obj_dict[:numlinks])]
  JuMP.set_optimizer(graph.obj_dict[:mflat], graph.optimizer)

  links = getlinkconstraints(graph)
  numlinks = length(links)

  for index in 1:num_nodes(graph)
    node = getnode(graph, index)
    model = getmodel(node)
    model.ext[:preobj] = objective_function(model)
    #Add theta to parent nodes
    if out_degree(graph,node) != 0
      childrenindices = [graph[child] for child in out_neighbors(graph,node)]
      sort!(childrenindices)
      @variable(model, θ[i in childrenindices] >= -1e6)
      set_objective_function(model, objective_function(model) + sum(θ[i] for i in childrenindices))
    end
  end
  #Add dual constraint to child nodes using the linking constraints
  for (numlink,link) in enumerate(links)
    #Take the two variables of the constraint
    var1, var2 = collect(keys(link.func.terms))
    #Determine which nodes they belong to
    nodeV1 = getnode(var1)
    nodeV2 = getnode(var2)
    #Set the order of the nodes
    if nodeV1 in out_neighbors(graph,nodeV2)
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
    childindex = graph[childnode]
    childmodel = getmodel(childnode)
    push!(parentnode[:childvars][childindex],parentvar)
    valbar = @variable(childmodel)
    set_name(valbar,"varlbar$numlink")
    push!(childnode[:xinvars],valbar)
    conref = @constraint(childmodel, valbar - childvar == 0)
    push!(childnode[:linkconstraints], conref)
  end
  graph.obj_dict[:preprocessed] = true
end

function sendnodestoworkers(graph)
  levels = graph.obj_dict[:levels]
  numlevels = length(levels)
  for level in 1:numlevels
    nodeslevel = levels[level]
    worker = repeat(workers(),outer=Int(ceil(length(nodeslevel)/nworkers())))
    for (k,node) in enumerate(nodeslevel)
      println("Sending node $(node.label) to worker $(worker[k])")  
      remotecall_fetch(()->node,worker[k])
    end
  end
 end    
include("liftproject.jl")
include("gomory.jl")
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
function bendersolve(graph::Plasmo.PlasmoGraph; max_iterations::Int64=10, cuts::Array{Symbol,1}=[:LP], ϵ=1e-5, rel_gap=1e-4, UBupdatefrequency=1,timelimit=3600,verbose=false, LB=NaN, is_nonconvex=false)
  graph.attributes[:is_nonconvex] = is_nonconvex
  global tmpdir = "/tmp/RootNode" # mktempdir()
  
  updatebound = true

  verbose && info("Preparing graph")
  bdprepare(graph, cuts)
  n = graph.attributes[:normalized]
  if isnan(LB) && (!is_nonconvex)
    verbose && info("Solve relaxation and set LB")
    mf = graph.attributes[:mflat]
    solve(mf,relaxation=true)
    LB = getobjectivevalue(graph.attributes[:mflat])
  end
  UB = Inf

  # Set bound to root node
  rootnode = graph.attributes[:roots][1]
  rootmodel = getmodel(rootnode)
  # @constraint(rootmodel, rootmodel.obj.aff >= LB)

  # Begin iterations
  verbose && info("Begin iterations")
  starttime = time()
  for i in 1:max_iterations
    tic()
    updatebound = ((i-1) % UBupdatefrequency) == 0
    LB,UB = forwardstep(graph, cuts, updatebound)

    tstamp = time() - starttime

    itertime = toc()
    if n == 1
      saveiteration(graph.attributes[:solution],tstamp,[graph.attributes[:iterUB],LB,itertime,tstamp],n)
    else
      saveiteration(graph.attributes[:solution],tstamp,[n*LB,n*UB,itertime,tstamp],n)
    end
    printiterationsummary(graph.attributes[:solution],singleline=false)

    if abs(UB-LB) < ϵ
      graph.attributes[:solution].termination = "Optimal"
      return graph.attributes[:solution]
    end

    if calculate_gap(LB, UB) < rel_gap
      graph.attributes[:solution].termination = "Optimal"
      return graph.attributes[:solution]
    end

    # Check time limit
    if tstamp > timelimit
      graph.attributes[:solution].termination = "Time Limit"
      return graph.attributes[:solution]
    end

    if graph.attributes[:stalled]
      graph.attributes[:solution].termination = "Stalled"
      return graph.attributes[:solution]
    end

    if graph.attributes[:fathomed_by_bound]
      graph.attributes[:solution].termination = "node fathomed by bound"
      return graph.attributes[:solution]
    end

    if graph.attributes[:is_infeasible]
      graph.attributes[:solution].termination = "Benders master problem is infeasible"
      return graph.attributes[:solution]
    end

  end

  graph.attributes[:solution].termination = "Max Iterations"
  return graph.attributes[:solution]
end

function forwardstep(graph::PlasmoGraph, cuts::Array{Symbol,1}, updatebound::Bool)
  levels = graph.attributes[:levels]
  numlevels = length(levels)
  for level in 1:numlevels
    nodeslevel = levels[level]
    for node in nodeslevel
      solveprimalnode(node,graph,cuts,updatebound)
    end
  end
  LB = graph.attributes[:LB]
  if updatebound
    iterUB = sum(node.attributes[:preobjval] for node in values(graph.nodes))
    graph.attributes[:iterUB] = iterUB
    UB = min(graph.attributes[:UB],iterUB)
    #update best feasible x
    if UB < graph.attributes[:UB]
      root = graph.attributes[:roots][1]
      if !haskey(graph.attributes, :best_feasible_x)
        graph.attributes[:best_feasible_x] = Dict()
      end
      for varname in keys(root.attributes[:varname_to_var])
        graph.attributes[:best_feasible_x][varname] = getvalue(root.attributes[:varname_to_var][varname])
      end
    end
    graph.attributes[:UB] = UB
  else
    UB = graph.attributes[:UB]
  end
  return LB,UB
end

function solveprimalnode(node::PlasmoNode, graph::PlasmoGraph, cuts::Array{Symbol,1}, updatebound::Bool)
  # 1. Add cuts
  generatecuts(node,graph)
  # 2. Take x
  takex(node, graph)
  # 3. solve
  if :LP in cuts && in_degree(graph,node) != 0
    solvelprelaxation(node)
  end

  if :GMI in cuts && in_degree(graph,node) != 0
    solvegmirelaxation(node, graph)
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

function solvelprelaxation(node::PlasmoNode)
  model = getmodel(node)
  status = solve(model, relaxation = true)

  # if status != :Optimal
  #   println(node.attributes[:xin])
  #   error("lp relaxation not solved to optimality")
  # end

  dualconstraints = node.attributes[:linkconstraints]

  λnode = getdual(dualconstraints)
  nodebound = getobjectivevalue(model)

  node.attributes[:bound] = nodebound
  node.attributes[:λ] = λnode

  return status
end




function solvenodemodel(node::PlasmoNode,graph::PlasmoGraph)
  if graph.attributes[:is_nonconvex] && in_degree(graph, node) != 0
    model = node.attributes[:ubsub]
    # println(model)
    status = solve(model)

    node.attributes[:preobjval] = getvalue(node.attributes[:preobj])
  else 
    model = getmodel(node)
    status = solve(model)
    if status != :Optimal && in_degree(graph, node) == 0
      graph.attributes[:is_infeasible] = true
      graph.attributes[:LB] = +Inf
    end
    if status != :Optimal && in_degree(graph, node) !=0
      println(node.attributes[:xin])
      # error("upper bound subproblem not solved to optimality")
    end
    if in_degree(graph,node) == 0 # Root node
      graph.attributes[:LB] = getobjectivevalue(model)
    end
    if graph.attributes[:LB] > graph.attributes[:global_UB]
      graph.attributes[:fathomed_by_bound] = true
    end
    node.attributes[:preobjval] = getvalue(node.attributes[:preobj])
  end

end

function takex(node::PlasmoNode, graph::PlasmoGraph)
  if in_degree(graph, node) == 0 
    return true
  end
  xinvals = node.attributes[:xin]

  xinvars = node.attributes[:xinvars]
  if length(xinvals) > 0
    fix.(xinvars,xinvals)
    if graph.attributes[:is_nonconvex]
      fix.(node.attributes[:ubxinvars], xinvals)
    end
  end
end

function putx(node::PlasmoNode,graph::PlasmoGraph)
  childvars = node.attributes[:childvars]
  children = out_neighbors(graph,node)
  length(children) == 0 && return true

  for child in children
    xnode = getvalue(childvars[getnodeindex(graph,child)])
    #round xnode to bounds (sometimes numerical errors can occur)
    for i in 1:length(childvars[getnodeindex(graph,child)])
      var = childvars[getnodeindex(graph,child)][i]
      ub = getupperbound(var)
      lb = getlowerbound(var)
      if xnode[i] > ub 
        xnode[i] = ub 
      end
      if xnode[i] < lb 
        xnode[i] = lb 
      end
    end
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
  if :LP in cuts || :GMI in cuts || :LIFT in cuts
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
  previouscuts = node.attributes[:prevcuts]
  thisitercuts = Dict()
  samecuts = Dict()
  for child in children
    childindex = getnodeindex(graph,child)
    thisitercuts[childindex] = CutData[]
    samecuts[childindex] = Bool[]

    while length(cutdataarray[childindex]) > 0
      cutdata = pop!(cutdataarray[childindex])
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

  cuts = graph.attributes[:cuts]
  if :LP in cuts 
    node.attributes[:stalled] = reduce(*,nodesamecuts)
    #bound does not improve for 5 iterations is also considered as stalled
    if length(graph.attributes[:solution].iterbound) > 6 && abs(graph.attributes[:LB] - graph.attributes[:solution].iterbound[end-5]) < 1e-3
      node.attributes[:stalled] = true
    end 
    node.attributes[:stalled] && warn("Node $(node.label) stalled")
  end

  #only stalled when LP has already stalled
  if (:GMI in cuts || :LIFT in cuts) && node.attributes[:LP_stalled]
    node.attributes[:stalled] = reduce(*,nodesamecuts)
    if length(graph.attributes[:solution].iterbound) > 6 && abs(graph.attributes[:LB] - graph.attributes[:solution].iterbound[end-5]) < 1e-3 && node.attributes[:LP_stalled_iterations] > 2
      node.attributes[:stalled] = true
    end    
    node.attributes[:LP_stalled_iterations] += 1 
    node.attributes[:stalled] && warn("Node $(node.label) stalled")
  end 
  if (:GMI in cuts || :LIFT in cuts) && (!node.attributes[:LP_stalled])
    node.attributes[:LP_stalled] = reduce(*,nodesamecuts)
    if length(graph.attributes[:solution].iterbound) > 6 && abs(graph.attributes[:LB] - graph.attributes[:solution].iterbound[end-5]) < 1e-3
      node.attributes[:LP_stalled] = true
    end      
    println("LP_stalled status")
    println(node.attributes[:LP_stalled])
    if node.attributes[:LP_stalled]
      node.attributes[:LP_stalled_iterations] = 1
    end 
  end 
  if in(node,graph.attributes[:roots]) && node.attributes[:stalled]
    graph.attributes[:stalled] = true
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
    node.attributes[:xin] = []
    node.attributes[:λ] = []
    node.attributes[:bound] = NaN
    node.attributes[:xinvars] = []
    node.attributes[:preobjval] = NaN
    node.attributes[:linkconstraints] = []
    if graph.attributes[:is_nonconvex] && in_degree(graph, node) !=0 
      node.attributes[:ubxinvars] = []
    end
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
        node.attributes[:prevcuts] = Dict(getnodeindex(graph,child) => CutData[] for child in out_neighbors(graph,node))
        node.attributes[:stalled] = false
    end
    current = children
    level += 1
  end
  graph.attributes[:numlevels] = level - 1
end

function bdprepare(graph::Plasmo.PlasmoGraph, cuts::Array{Symbol,1}=[:LP])
  if haskey(graph.attributes,:preprocessed)
    return true
  end
  graph.attributes[:solution]= Solution(method=:benders)
  identifylevels(graph)
  graph.attributes[:normalized] = normalizegraph(graph)
  graph.attributes[:stalled] = false
  if !graph.attributes[:is_nonconvex]
    graph.attributes[:mflat] = create_flat_graph_model(graph)
    setsolver(graph.attributes[:mflat],graph.solver)
  end
  graph.attributes[:UB] = Inf
  graph.attributes[:global_UB] = +Inf
  graph.attributes[:fathomed_by_bound] = false
  graph.attributes[:is_infeasible] = false
  graph.attributes[:cuts] = cuts

  links = getlinkconstraints(graph)
  numlinks = length(links)

  for index in 1:length(graph.nodes)
    node = graph.nodes[index]
    model = getmodel(node)
    if model.solver == JuMP.UnsetSolver()
      model.solver = graph.solver
    end
    if graph.attributes[:is_nonconvex] && in_degree(graph, node) != 0
      node.attributes[:preobj] = node.attributes[:ubsub].obj 
    else
      node.attributes[:preobj] = model.obj
    end

    
    if :GMI in cuts || :LIFT in cuts 
    #create upper bound for variables if GMI or LIFT is in cuts   
      for col in 1:length(model.colUpper)
        if model.colUpper[col] > 1e10 
          setupperbound(Variable(model, col), 1e10)
        end
      end
      #add attributes to check if LP cuts has stalled 
      node.attributes[:LP_stalled] = false 
    end

    #create standard matrix for lift and project cuts 
    if :LIFT in cuts && in_degree(graph, node) != 0
      getmatrixform(node)
    end

    #map variable name to variable 
    node.attributes[:varname_to_var] = Dict()
    node.attributes[:varname_to_varcat] = Dict()
    for i in 1:length(model.colNames)
      varname = model.colNames[i]
      node.attributes[:varname_to_var][varname] = Variable(model, i)
      node.attributes[:varname_to_varcat][varname] = model.colCat[i]
    end

    if graph.attributes[:is_nonconvex] && in_degree(graph, node) != 0
      node.attributes[:nonconvex_varname_to_var] = Dict()
      for i in 1:length(node.attributes[:ubsub].colNames)
        varname = node.attributes[:ubsub].colNames[i]
        node.attributes[:nonconvex_varname_to_var][varname] = Variable(node.attributes[:ubsub], i)
      end
    end

    if out_degree(graph,node) != 0
#Add theta to parent nodes
      childrenindices = [getnodeindex(graph,child) for child in out_neighbors(graph,node)]
      sort!(childrenindices)
      @variable(model, θ[i in childrenindices] >= -1e6)
      model.obj += sum(θ[i] for i in childrenindices)
      node.attributes[:θ_to_scenario] = Dict()
      node.attributes[:scenario_to_θ] = Dict()
      for child in out_neighbors(graph,node)
        if haskey(child.attributes, :scenario)
          index = getnodeindex(graph, child)
          node.attributes[:θ_to_scenario][index] = child.attributes[:scenario]
          node.attributes[:scenario_to_θ][child.attributes[:scenario]] = index
        end
      end

    end

    #get the index of linking variables 
    if in_degree(graph, node) != 0
      node.attributes[:linking_vars_indices] = []
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
    #store child var index 
    push!(childnode.attributes[:linking_vars_indices], childvar.col)
    push!(childnode.attributes[:linking_vars_indices], valbar.col)
  end

  #link master and ub subproblem for nonconvex 
  if graph.attributes[:is_nonconvex]
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
      childmodel = childnode.attributes[:ubsub]
      # println(childmodel)
      # println(childmodel.colNames)
      temp_model = childvar.m 
      childvar_name = temp_model.colNames[childvar.col]
      # println(childvar_name)
      # println(get_col_from_varname(childmodel, childvar_name))
      # push!(parentnode.attributes[:childvars][childindex],parentvar)
      valbar = @variable(childmodel)
      setname(valbar,"varlbar$numlink")
      push!(childnode.attributes[:ubxinvars],valbar)
      @constraint(childmodel, valbar - Variable(childmodel,get_col_from_varname(childmodel, childvar_name)) == 0)
    end
  end

  #set up CGLPs
  if :LIFT in cuts 
    for index in 1:length(graph.nodes)
      node = graph.nodes[index]
      if in_degree(graph, node) != 0 
        setCGLP(node, graph)
      end
    end
  end

  graph.attributes[:preprocessed] = true
end



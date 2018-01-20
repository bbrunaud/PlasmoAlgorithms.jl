using JuMP
using Gurobi
using Plasmo

function fix(var,value)
  ##Sets value for constraint variable
  setlowerbound(var,value)
  setupperbound(var,value)
end

function isChildNode(g::Plasmo.PlasmoGraph, n1::PlasmoNode, n2::PlasmoNode)
  ##Checks if n1 is a child node of n2
  for node in LightGraphs.out_neighbors(g.graph,getindex(g,n2))
    if (n1 == g.nodes[node]) return true
      return true
    end
  end
  return false
end

function numChildNodes(g::PlasmoGraph, n1::PlasmoNode)
  return length(LightGraphs.out_neighbors(g.graph,getindex(g,n1)))
end

function numParentNodes(g::PlasmoGraph, n1::PlasmoNode)
  return length(LightGraphs.in_neighbors(g.graph,getindex(g,n1)))
end

function levels(graph::Plasmo.PlasmoGraph)
  #Create lists of root and leaf nodes in graph
  graph.attributes[:roots] = []
  graph.attributes[:leaves] = []
  #Create dictionary to keep track of levels of nodes
  graph.attributes[:levels] = Dict()
  #Iterate through every node to check for root/leaf nodes
  for nodeIndex in 1:length(graph.nodes)
    node = graph.nodes[nodeIndex]
    #If the node does not have parents it is a root node
    if numParentNodes(graph,node)==0
      push!(graph.attributes[:roots],node)
      #Root nodes are the first level
    end
    #If the node does not have children it is a leaf node
    if numChildNodes(graph,node)==0
      push!(graph.attributes[:leaves],node)
    end
  end

  #Start mapping level from the root nodes
  currentNodes = graph.attributes[:roots]
  levels = graph.attributes[:levels]
  level = 1
  while currentNodes != []
    levels[level] = currentNodes
    childrenNodes = []
    for node in currentNodes
      childNodeIndex = LightGraphs.out_neighbors(graph.graph,getindex(graph,node))
      for index in childNodeIndex
        childNode = graph.nodes[index]
        push!(childrenNodes,childNode)
      end
    end
    level += 1
    currentNodes = childrenNodes
  end
end

function preProcess(graph::Plasmo.PlasmoGraph)
  levels(graph)

  links = getlinkconstraints(graph)
  numLinks = length(links)
  numNodes = length(graph.nodes)

  graph.attributes[:links] = Dict()
  graph.attributes[:childlinks] = Dict()
  graph.attributes[:duals] = Dict()
  dualMap = graph.attributes[:duals]
  linkIndex = graph.attributes[:links]
  childLinkIndex = graph.attributes[:childlinks]

  #Add each node as a key to the dictionary
  for index in 1:length(graph.nodes)
    node = graph.nodes[index]
    nodeModel = getmodel(node)
    nodeModel.ext[:preobj] = nodeModel.obj
    linkIndex[node] = []
    childLinkIndex[node] = []
    #Add valbar to child nodes
    if numParentNodes(graph,node) != 0
      sp = getmodel(node)
      @constraintref dual[1:numLinks]
      @variable(sp, valbar[1:numLinks])
    end
    #Add theta to parent nodes
    if numChildNodes(graph,node) != 0
      mp = getmodel(node)
      mflat = create_flat_graph_model(graph)
      bound = getobjectivevalue(mflat)
      @variable(mp,θ[1:numNodes] >= -1e6)
      for node in LightGraphs.out_neighbors(graph.graph,getindex(graph,node))
        childNode = graph.nodes[node]
        childNodeIndex = getindex(graph,childNode)
        θ = getindex(mp,:θ)
        mp.obj += θ[childNodeIndex]
      end
    end
  end

  #Add dual constraint to child nodes using the linking constraints
  for link in 1:numLinks
    #Take the two variables of the constraint
    var1 = links[link].terms.vars[1]
    var2 = links[link].terms.vars[2]
    #Determine which nodes they belong to
    nodeV1 = getnode(var1)
    nodeV2 = getnode(var2)
    #Set the order of the nodes
    if isChildNode(graph,nodeV1,nodeV2)
      childNode = nodeV1
      childvar = var1
      parentNode = nodeV2
    elseif isChildNode(graph,nodeV2,nodeV1)
      childNode = nodeV2
      childvar = var2
      parentNode = nodeV1
    end
    #Add theta[childNode] to the parent node objective
    mp = getmodel(parentNode)
    childNodeIndex = getindex(graph,childNode)
    # θ = getindex(mp,:θ)
    # mp.obj += θ[childNodeIndex]
    #Add linking constraint to the parent node dictionary
    linkList = linkIndex[parentNode]
    childLinkList = childLinkIndex[childNode]
    push!(linkList,link)
    push!(childLinkList,link)
  end
  for child in keys(childLinkIndex)
    sp = getmodel(child)
    childLinkList = childLinkIndex[child]
    @constraintref dual[1:length(childLinkList)]
    for i in 1:length(childLinkList)
      link = childLinkList[i]
      var1 = links[link].terms.vars[1]
      var2 = links[link].terms.vars[2]
      #Determine which nodes they belong to
      nodeV1 = getnode(var1)
      nodeV2 = getnode(var2)
      #Set the order of the nodes
      if isChildNode(graph,nodeV1,nodeV2)
        childNode = nodeV1
        childvar = var1
        parentNode = nodeV2
      elseif isChildNode(graph,nodeV2,nodeV1)
        childNode = nodeV2
        childvar = var2
        parentNode = nodeV1
      end
      valbar = getindex(sp,:valbar)
      dual[i] = @constraint(sp,valbar[link] - childvar == 0)
    end
    dualMap[child] = dual
  end
end

function forwardStep(graph::Plasmo.PlasmoGraph)
  levels = graph.attributes[:levels]
  links = getlinkconstraints(graph)
  linksMap = graph.attributes[:links]
  for level in 1:length(levels)
    currentLevel = levels[level]
    for node in currentLevel
      mp = getmodel(node)
      print(mp)
      solve(mp)
      #Get the constraints linked to this node from the dictionary
      nodelinks = linksMap[node]
      #Iterate through linking constraints
      for link in nodelinks
        #Get the nodes and variables in the linked constraint
        var1 = links[link].terms.vars[1]
        var2 = links[link].terms.vars[2]
        nodeV1 = getnode(var1)
        nodeV2 = getnode(var2)
        #Determine which nodes are parents and children
        if isChildNode(graph,nodeV1,nodeV2)
          childNode = nodeV1
          parentNode = nodeV2
          val = getvalue(var2)
        elseif isChildNode(graph,nodeV2,nodeV1)
          childNode = nodeV2
          parentNode = nodeV1
          val = getvalue(var1)
        end
        #Set valbar in child node associated with link to variable value in parent
        sp = getmodel(childNode)
        valbar = getindex(sp, :valbar)
        fix(valbar[link], val)
      end
    end
  end

  leaves = graph.attributes[:leaves]
  roots = graph.attributes[:roots]
  LB = 0
  UB = 0
  for root in roots
    rootmodel = getmodel(root)
    LB = LB + getobjectivevalue(rootmodel)
    # UB += getvalue(rootmodel.ext[:preobj])
  end
  for index in 1:length(graph.nodes)
    node = graph.nodes[index]
    nodeModel = getmodel(node)
    UB += getvalue(nodeModel.ext[:preobj])
  end
  return LB,UB
end

function backwardStep(cut::String,graph::Plasmo.PlasmoGraph)
  levels = graph.attributes[:levels]
  for level in length(keys(levels)):-1:2
    for node in levels[level]
      cutGeneration(cut,graph,node)
    end
  end
end

function bendersolve(cut::String, graph::Plasmo.PlasmoGraph; max_iterations = 3)
  preProcess(graph)
  ϵ = 10e-5
  UB = Inf
  LB = -Inf

  for i in 1:max_iterations
    LB,UB = forwardStep(graph)
    println("*** ",UB)
    println("*** ",LB)
    # run(`cowsay "UB = $UB"`)
    if abs(UB - LB)<ϵ
      println("Converged!")
      # run(`cowsay -f moofasa "Converged!"`)
      break
    end
    backwardStep(cut,graph)
  end
  return getobjectivevalue(getmodel(graph.nodes[1]))
end

function cutGeneration(cut::String, graph::PlasmoGraph, node::PlasmoNode)
  if cut == "LP"
    LPcut(graph,node)
  end
end

function LPcut(graph::PlasmoGraph, node::PlasmoNode)
  linkList= graph.attributes[:links]
  dualMap = graph.attributes[:duals]
  childLinks = graph.attributes[:childlinks]
  links = getlinkconstraints(graph)

  sp = getmodel(node)
  λs = []
  valbars = []
  variables = []
  solve(sp)
  for childLink in childLinks[node]
    # dualCon = getindex(sp,:dual)
    dualCon = dualMap[node]
    λs = getdual(dualCon)
    valbar = getindex(sp,:valbar)
    println(sp)
    push!(valbars,valbar[childLink])

    var1 = links[childLink].terms.vars[1]
    var2 = links[childLink].terms.vars[2]
    nodeV1 = getnode(var1)
    nodeV2 = getnode(var2)
    #Determine which nodes are parents and children
    if isChildNode(graph,nodeV1,nodeV2)
      var = var2
    elseif isChildNode(graph,nodeV2,nodeV1)
      var = var1
    end
    push!(variables,var)
  end
  parentNodes = LightGraphs.in_neighbors(graph.graph,getindex(graph,node))
  parentNode = graph.nodes[parentNodes[1]]

  mp = getmodel(parentNode)
  θ = getindex(mp,:θ)

  status = solve(getmodel(node), relaxation = true)
  rhs = 0

  for i in 1:length(variables)
    rhs = rhs + λs[i]*(getupperbound(valbars[i])-variables[i])
  end
  if status != :Optimal
    @constraint(mp, 0 >= rhs)
    println("HERE")
  else
    θk = getobjectivevalue(getmodel(node))
    @constraint(mp, θ[getindex(graph,node)] >= θk + rhs)
  end
end

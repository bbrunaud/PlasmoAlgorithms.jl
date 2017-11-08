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

function levelsRecursive(g::Plasmo.PlasmoGraph,n::PlasmoNode,level::Int64,levelDict::Dict)
  if level in keys(levelDict)
  else
    levelDict[level] = Set()
  end
  levelInv = levelDict[level]
  push!(levelInv,n)
  if numChildNodes(g,n) == 0
    #lastLevel
    return 1
  end
  for index in LightGraphs.out_neighbors(g.graph,getindex(g,n))
    childNode = g.nodes[index]
    levelsRecursive(g,childNode,level+1,levelDict)
  end
end

function getLevels(g::Plasmo.PlasmoGraph)
  levelDict = Dict()
  levelCount = 1
  currentNode = g.nodes[1]
  levelsRecursive(g,currentNode,levelCount, levelDict)
  return levelDict
end

function preProcess(graph::Plasmo.PlasmoGraph)
  #Build a dictionary for the different levels of the tree and each node in it
  levels = getLevels(graph)

  #Grab the linking constraints from the graph
  links = getlinkconstraints(graph)
  numLinks = length(links)
  numNodes = length(graph.nodes)

  #Create a dictionary for node relationships
  dict = Dict()
  childDict = Dict()

  #Add each node as a key to the dictionary
  for index in 1:length(graph.nodes)
    node = graph.nodes[index]
    dict[node] = Set()
    #Add valbar to child nodes
    if numParentNodes(graph,node) != 0
      sp = getmodel(node)
      @variable(sp, valbar[1:numLinks])
    end
    #Add theta to parent nodes
    if numChildNodes(graph,node) != 0
      mp = getmodel(node)
      @variable(mp,θ[1:numNodes] >= 0)
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
    θ = getindex(mp,:θ)
    mp.obj += θ[childNodeIndex]
    #Add dual constraint to the child node model
    sp = getmodel(childNode)
    valbar = getindex(sp,:valbar)
    parentNodeIndex = getindex(graph,parentNode)
    @constraint(sp, dual, valbar[link] - childvar == 0)
    print(sp)
    #Add linking constraint to the parent node dictionary
    linkIndex = dict[parentNode]
    push!(linkIndex,link)
  end
  return dict, levels
end

function forwardStep(graph::Plasmo.PlasmoGraph, dict::Dict, node::Plasmo.PlasmoNode)
  #If node has no children exit function; we have reached the end of the tree
  if numChildNodes(graph, node) == 0
    return 0
  end
  #Solve the current node
  mp = getmodel(node)
  solve(mp)
  #Get the constraints linked to this node from the dictionary
  links = getlinkconstraints(graph)
  nodelinks = dict[node]
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
    parentIndex = getindex(g,parentNode)
    valbar = getindex(sp, :valbar)
    fix(valbar[link], val)
    #Repeat this step for the child node
    forwardStep(graph, dict, childNode)
  end
  for index in 1:length(graph.nodes)
    testNode = graph.nodes[index]
    model = getmodel(testNode)
  end
end

function backwardStep(graph::Plasmo.PlasmoGraph, levels::Dict, dict::Dict)
  for key in length(keys(levels)):-1:1
    for node in levels[key]
      if numParentNodes(graph,node) != 0
        parentNodes = LightGraphs.in_neighbors(graph.graph,getindex(graph,node))
        parentNode = graph.nodes[parentNodes[1]]
        parentIndex = getindex(graph, parentNode)
        childIndex = getindex(graph, node)
        sp = getmodel(node)
        mp = getmodel(parentNode)
        status = solve(sp)
        dualCon = getindex(sp,:dual)
        λ = getdual(dualCon)
        debug("Dual solution is found to be $λ")

        links = getlinkconstraints(graph)
        nodeLinks = dict[parentNode]
        for link in nodeLinks
          #Get the nodes and variables in the linked constraint
          var1 = links[link].terms.vars[1]
          var2 = links[link].terms.vars[2]
          nodeV1 = getnode(var1)
          nodeV2 = getnode(var2)
          #Determine which nodes are parents and children
          if isChildNode(graph,nodeV1,nodeV2)
            childNode = nodeV1
            parentNode = nodeV2
            var = var2
            val = getvalue(var2)
          elseif isChildNode(graph,nodeV2,nodeV1)
            childNode = nodeV2
            parentNode = nodeV1
            val = getvalue(var1)
            var = var1
          end

          if node == childNode
            valbar = getindex(sp,:valbar)
            θ = getindex(mp,:θ)
            if status != :Optimal
              @constraint(mp, 0 >= λ*(getupperbound(valbar[parentIndex])-var))
              println(mp)
            else
              θk = getobjectivevalue(sp)
              @constraint(mp, θ[childIndex] >= θk + λ*(getvalue(valbar[parentIndex])-var))
            end
          end
        end
      end
    end
  end
  print(mp)
end

function bendersolve(graph::Plasmo.PlasmoGraph, max_iterations::Int64)
    dict, levels = preProcess(graph)
    numLevels = length(levels)
    for i in 1:max_iterations
      debug("Iteration $i")
      forwardStep(graph, dict, graph.nodes[1])
      LB = getobjectivevalue(getmodel(graph.nodes[1]))
      debug("Lowerbound = $LB")
      ##TODO ask Braulio about where the UB needs to be
      UB = 0
      if UB == LB
        break
      end
      backwardStep(graph, levels, dict)
    end
    return getobjectivevalue(getmodel(graph.nodes[1]))
end

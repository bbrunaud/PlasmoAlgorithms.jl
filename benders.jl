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

function bendersetup(graph::PlasmoGraph)
  ##Add all linked constraint to child node models
  links = getlinkconstraints(graph)
  numLinks = length(links)
  dict = Dict()

  for node in 1:length(graph.nodes)
    dict[graph.nodes[node]] = Set()
  end

  for link in 1:numLinks
    var1 = links[link].terms.vars[1]
    var2 = links[link].terms.vars[2]

    nodeV1 = getnode(var1)
    nodeV2 = getnode(var2)
    if isChildNode(graph,nodeV1,nodeV2)
      sp = getmodel(nodeV1)
      linkIndex = dict[nodeV2]
      push!(linkIndex,link)
      @variable(sp, valbar[1:numLinks])
    elseif isChildNode(graph,nodeV2,nodeV1)
      sp = getmodel(nodeV2)
      linkIndex = dict[nodeV1]
      push!(linkIndex,link)
      @variable(sp, valbar[1:numLinks])
    end
  end
  return [graph,dict]
end

function benderparent(graph::PlasmoGraph, dict::Dict, max_iterations::Int64, currentNode::PlasmoNode,nodeIndex)
  numParents = numParentNodes(graph, currentNode)
  if numParents == 0
    return getobjectivevalue(getmodel(currentNode))
  end
  
  parentNodes = LightGraphs.in_neighbors(g.graph,getindex(g,currentNode))
  parentNodeIndex = parentNodes[1]
  parentNode = graph.nodes[parentNodeIndex]
  
  numNodes = length(graph.nodes)
  
  mp = getmodel(parentNode)
  sp = getmodel(currentNode)

  #TODO change to flattening graph and adding bound
  @variable(mp,θ[1:numNodes])
  @constraint(mp,θ[nodeIndex]>=0)
  mp.obj += θ[nodeIndex]

  for i in max_iterations
    solve(mp)

    links = getlinkconstraints(graph)
    nodelinks = dict[parentNode]

    for link in nodelinks
      var1 = links[link].terms.vars[1]
      var2 = links[link].terms.vars[2]

      nodeV1 = getnode(var1)
      nodeV2 = getnode(var2)

      if isChildNode(graph,nodeV1,nodeV2)
        sp = getmodel(nodeV1)
        val = getvalue(var2)
        fix(valbar[link],val)
        @constraint(sp, dual, valbar[link] - var1 == 0)
      elseif isChildNode(graph,nodeV2,nodeV1)
        sp = getmodel(nodeV2)
        val = getvalue(var1)
        println("*help")
        #TODO ask Braulio why this doesn't recognize valbar as a variable
        fix(valbar[link],val)
        @constraint(sp, dual, valbar[link] - var2 == 0)
      end
      
      status = solve(sp)
      λ = getdual(dual)

      if status != :Optimal
        @constraint(mp, 0>=λ*(getupperbound(valbar)-var))
        println(mp)
      else
        θk = getobjectivevalue(sp)
        if θk == getvalue(θ[nodeIndex])
          benderparent(graph,dict,max_iterations, parentNode)
        end
        @constraint(mp, θ[nodeIndex] >= θk + λ*(getvalue(valbar)-var))
      end
    end
  end
end

function bendersrecursive(graph::PlasmoGraph, dict::Dict, max_iterations::Int64, nodeIndex = 1)
  currentNode = graph.nodes[nodeIndex]
  numChildren = numChildNodes(graph, currentNode)
  if numChildren == 0
    #do benders
    benderparent(graph,dict,max_iterations,currentNode,nodeIndex)
  else
    childrenIndex = LightGraphs.out_neighbors(g.graph,getindex(g,currentNode))
    for child in 1:length(childrenIndex)
      bendersrecursive(graph,dict,max_iterations,childrenIndex[child])
    end
  end
end

function bendersolve(graph::Plasmo.PlasmoGraph; max_iterations = 10)
    prepWork = bendersetup(graph)
    g = prepWork[1]
    dict = prepWork[2]
    bendersrecursive(g,dict,max_iterations)
end

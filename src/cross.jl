using JuMP
using Plasmo

function getParentNode(graph::PlasmoGraph, node::PlasmoNode)
  return LightGraphs.in_neighbors(graph.graph,getindex(graph, node))
end

function numParentNodes(graph::PlasmoGraph, node::PlasmoNode)
  return length(LightGraphs.in_neighbors(graph.graph,getindex(graph, node)))
end

function auxGraph(graph::PlasmoGraph, node::PlasmoNode)
  aGraph = PlasmoGraph()
  aGraph.solver = graph.solver
  add_node!(aGraph, node)
  origNodes = [getindex(graph,node)]
  while numParentNodes(graph, node) > 0
    parents = getParentNode(graph, node)
    for parent in parents
        add_node!(aGraph, graph.nodes[parent])
        add_edge(aGraph, graph.nodes[parent], node)
        node = graph.nodes[parent]
        push!(origNodes, parent)
    end
  end
  flatModel = create_flat_graph_model(aGraph)
  return flatModel, origNodes
end

function addLinks(graph::PlasmoGraph)
  for nodeIndex in 1:(length(graph.nodes)-1)
    m1 = getmodel(graph.nodes[nodeIndex])
    for variable in m1.colNames
      for comparisonIndex in (nodeIndex + 1):length(graph.nodes)
        m2 = getmodel(graph.nodes[comparisonIndex])
        if variable in m2.colNames
          n1 = graph.nodes[nodeIndex]
          n2 = graph.nodes[comparisonIndex]
          @linkconstraint(graph,n1[Symbol(variable)]==n2[Symbol(variable)])
        end
      end
    end
  end
end

function crossPrepare(graph::PlasmoGraph)
  lGraph = PlasmoGraph()
  lGraph.solver = graph.solver
  lGraph.attributes[:bGraph] = graph

  identifylevels(graph)
  leaves = graph.attributes[:leaves]

  for leaf in leaves
    flatModel, origNodes = auxGraph(graph, leaf)
    n = add_node(lGraph)
    setmodel(n,flatModel)
  end
  addLinks(lGraph)
  bdprepare(lGraph.attributes[:bGraph])
  return lGraph
end

function crossSolve(graph::PlasmoGraph)
  lGraph = crossPrepare(graph)
  lagrangesolve(lGraph, max_iterations = 1, update_method=:cuttingplanes, lagrangeheuristic = crossheuristic)
end

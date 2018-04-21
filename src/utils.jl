function normalizegraph(graph::PlasmoGraph)
    n = 1
    nodedict = getnodes(graph)
    for k in keys(nodedict)
        node = nodedict[k]
        m = getmodel(node)
        if m.objSense == :Max
            m.objSense = :Min
            m.obj = -m.obj
            n = -1
        end
    end
    graph.attributes[:normalized] = n
    return n
end

function fix(var::JuMP.Variable,value::Real)
  setlowerbound(var,value)
  setupperbound(var,value)
end

"""
Checks if n1 is a child node of n2
"""
ischildnode(graph::PlasmoGraph, n1::PlasmoNode, n2::PlasmoNode) = in(n2,in_neighbors(graph,n1))
 

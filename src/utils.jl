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

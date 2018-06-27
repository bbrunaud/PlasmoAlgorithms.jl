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

function savenodeobjective(mf::JuMP.Model)
    g = mf.ext[:Graph]
    numnodes = length(getnodes(g))
    nodeindex = Dict("node$i" => i for i in 1:numnodes)
    nov = mf.ext[:nodeobj] = [AffExpr(0.0) for i in 1:numnodes]
    obj = mf.obj.aff
    for (i,var) in enumerate(obj.vars)
        coeff = obj.coeffs[i]
        varname = mf.colNames[var.col]
        nodename = varname[1:search(varname,'.')-1]
        index = nodeindex[nodename]
        push!(nov[index],coeff,var)
    end
end

function get_frac(a)
    if length(a) == 1
        if abs(a%1) < 1e-10 || abs(a%1) > 1-1e-10
            return 0.0
        end
        if a > 0
            return a%1
        else
            return 1+(a%1)
        end
    end

    for i in 1:length(a)
        if abs(a[i]%1)<1e-10 || abs(a[i]%1) > 1-1e-10
            a[i] = 0.0
        else
            if a[i] > 0
                a[i] = a[i]%1
            else
                a[i]= 1+(a[i]%1)
            end
        end 
    end
    return a
end
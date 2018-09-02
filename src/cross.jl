mutable struct CrossGraph
 bd
 lg
end

function CrossGraph(g::ModelGraph)
 c = CrossGraph(g,deepcopy(g))
 return c
end

function crossprepare(c::CrossGraph)
    # cross-specific mappings
    # Lagrange mappings in Benders Graph
    links = getlinkconstraints(c.bd)
    for n in getnodes(c.bd)
        mn = getmodel(n)
        mn.ext[:preobj] = mn.obj
        mn.ext[:multmap] = Dict()
        mn.ext[:varmap] = Dict()
    end
    for (i,lc) in enumerate(links)
        for j in 1:length(lc.terms.vars)
            var = lc.terms.vars[j]
            var.m.ext[:multmap][i] = (lc.terms.coeffs[j],lc.terms.vars[j])
            var.m.ext[:varmap][var] = (i,j)
        end
    end

    # Benders mappings in Lagrange Graph
    PlasmoAlgorithms.identifylevels(c.lg)
    links = getlinkconstraints(c.lg)
    for (numlink,link) in enumerate(links)
        #Take the two variables of the constraint
        var1 = link.terms.vars[1]
        var2 = link.terms.vars[2]
        #Determine which nodes they belong to
        nodeV1 = getnode(var1)
        nodeV2 = getnode(var2)
        #Set the order of the nodes
        if ischildnode(c.lg,nodeV1,nodeV2)
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
        childindex = getindex(c.lg,childnode)
        push!(getattribute(parentnode, :childvars)[childindex],parentvar)
        if !hasattribute(parentnode, :λcomponents)
            setattribute(parentnode, :λcomponents, Dict())
        end
        if !haskey(getattribute(parentnode, :λcomponents), childindex)
            getattribute(parentnode, :λcomponents)[childindex] = Int64[]
        end
        push!(getattribute(parentnode, :λcomponents)[childindex],numlink)
    end
  
    # Standard preparations
    PlasmoAlgorithms.bdprepare(c.bd)
    PlasmoAlgorithms.lgprepare(c.lg, 0.5, 3, 1e6)
 end

function crosssolve(c::CrossGraph; max_iterations=5)
    crossprepare(c)
    for i in 1:max_iterations
        bendersolve(c.bd, max_iterations=1)
        benders_to_lagrange(c)
        lagrangesolve(c.lg,max_iterations=1, update_method=:cuttingplanes, callback=lagrange_to_benders)
        lagrange_to_benders(c)
    end
end

function benders_to_lagrange(c::CrossGraph)
    # Push Cut Data
    cutdataarray = getattribute(c.lg,:cutdata)
    for node in getnodes(c.bd)
        nodeindex = getnodeindex(node)
        m = getmodel(node)
        preobjval = JuMP.getvalue(m.ext[:preobj])
        λcomponent = Int64[]
        coeffs = Float64[]
        xk = Float64[]
        for k in keys(m.ext[:multmap])
            coeff = m.ext[:multmap][k][1]
            var = m.ext[:multmap][k][2]
            push!(λcomponent, k)
            push!(coeffs, coeff)
            push!(xk, JuMP.getvalue(var))
        end  
        push!(cutdataarray[nodeindex], LagrangeCutData(preobjval, λcomponent, coeffs, xk))
    end
end

function lagrange_to_benders(c::CrossGraph)
    λ = getattribute(c.lg,:λ)[end]
    for lgnode in getnodes(c.lg)
        if in_degree(c.lg, lgnode) > 0
            parentnode = in_neighbors(c.lg, lgnode)[1]
            nodeindex = PlasmoAlgorithms.getnodeindex(lgnode)
            parentindex = PlasmoAlgorithms.getnodeindex(parentnode)
            λk = λ[getattribute(parentnode,:λcomponents)[nodeindex]]
            bdnode = getnode(c.bd, parentindex)
            zld = PlasmoAlgorithms.subgraphobjective(lgnode, c.lg)
            push!(getattribute(bdnode,:cutdata)[nodeindex], LagrangeCrossCutData(zld,λk))
        end
    end
end

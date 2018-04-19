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

getx(graph::PlasmoGraph) = graph.attributes[:x]
getλ(graph::PlasmoGraph) = graph.attributes[:λ]
getrootnode(graph::PlasmoGraph) = graph.attributes[:rootnode]
getnodebounds(graph::PlasmoGraph) = graph.attributes[:nodebounds]
getnodeobjectivevals(graph::PlasmoGraph) = graph.attributes[:nodeobjectivevals]
getlevels(graph::PlasmoGraph) = graph.attributes[:levels]

function turnoffpresolve(model::JuMP.Model)
    if typeof(model.solver) == CPLEX.CplexSolver
        push!(model.solver.options,(:CPX_PARAM_PREIND,0))
    elseif typeof(model.solver) == Gurobi.GurobiSolver
        push!(model.solver.options,(:Presolve,0))
    else
        error("Root relaxation requires to turn off presolve for the solver, only CPLEX and Gurobi are supported")
    end
end

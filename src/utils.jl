
setattribute(graph::OptiGraph, attribute, value) = graph.obj_dict[attribute] = value
getattribute(graph::OptiGraph, attribute) = graph.obj_dict[attribute]
hasattribute(graph::OptiGraph, attribute) = haskey(graph.obj_dict, attribute)

setattribute(node::OptiNode, attribute, value) = node.model.obj_dict[attribute] = value
getattribute(node::OptiNode, attribute) = node.model.obj_dict[attribute]
hasattribute(node::OptiNode, attribute) = haskey(node.model.obj_dict, attribute)

function out_neighbors(graph::OptiGraph, node::OptiNode)
    neighbors = []
    edges = incident_edges(graph, [node])
    for e in edges
        e.nodes[1] == node && push!(neighbors, e.nodes[2])
    end

    return neighbors
end
out_degree(graph::OptiGraph, node::OptiNode) = length(out_neighbors(graph, node))


function in_neighbors(graph::OptiGraph, node::OptiNode)
    neighbors = []
    edges = incident_edges(graph, [node])
    for e in edges
        e.nodes[2] == node && push!(neighbors, e.nodes[1])
    end

    return neighbors
end
in_degree(graph::OptiGraph, node::OptiNode) = length(in_neighbors(graph, node))

function normalizegraph(graph::OptiGraph)
    n = 1
    for node in getnodes(graph)
        m = getmodel(node)
        if objective_sense(m) == MOI.MAX_SENSE
            set_objective_sense(m, MOI.MIN_SENSE)
            set_objective_function(m, -objective_function(m))
            n = -1
        end
    end
    setattribute(graph, :normalized, n)
    return n
end

function subgraphobjective(node, graph)
    if out_degree(graph, node) == 0
        return objective_value(getmodel(node))
    end
    return objective_value(getmodel(node)) + sum(subgraphobjective(childnode, graph) for childnode in out_neighbors(graph, node))
end

function relax_integrality(model::Model)
    info_pre_relaxation = map(v -> (v, _info_from_variable(v)),
        all_variables(model))
    # We gather the info first because some solvers perform poorly when you
    # interleave queries and changes. See, e.g.,
    # https://github.com/jump-dev/Gurobi.jl/pull/301.
    for (v, info) in info_pre_relaxation
        if info.integer
            unset_integer(v)
        elseif info.binary
            unset_binary(v)
            set_lower_bound(v, max(0.0, info.lower_bound))
            set_upper_bound(v, min(1.0, info.upper_bound))
        end
    end
    function unrelax()
        for (v, info) in info_pre_relaxation
            if info.integer
                set_integer(v)
            elseif info.binary
                set_binary(v)
                if info.has_lb
                    set_lower_bound(v, info.lower_bound)
                else
                    delete_lower_bound(v)
                end
                if info.has_ub
                    set_upper_bound(v, info.upper_bound)
                else
                    delete_upper_bound(v)
                end
            end
        end
        return
    end
    return unrelax
end

function _info_from_variable(v::VariableRef)
    has_lb = has_lower_bound(v)
    lb = has_lb ? lower_bound(v) : -Inf
    has_ub = has_upper_bound(v)
    ub = has_ub ? upper_bound(v) : Inf
    has_fix = is_fixed(v)
    fixed_value = has_fix ? fix_value(v) : NaN
    start_or_nothing = start_value(v)
    has_start = !(start_or_nothing isa Nothing)
    start = has_start ? start_or_nothing : NaN
    has_start = start !== Nothing
    binary = is_binary(v)
    integer = is_integer(v)
    return VariableInfo(has_lb, lb, has_ub, ub, has_fix, fixed_value,
                        has_start, start, binary, integer)
end 

function update!(rr::RemoteChannel,value)
    !isready(rr) && error("Updating an empty channel")
    take!(rr)
    put!(rr,value)
end

#=
import Plasmo.getgraph



function fix(var::JuMP.Variable,value::Real)
  setlowerbound(var,value)
  setupperbound(var,value)
end

"""
Checks if n1 is a child node of n2
"""
ischildnode(graph::ModelGraph, n1::ModelNode, n2::ModelNode) = in(n2,in_neighbors(graph,n1))

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

function getnodeindex(node::Plasmo.PlasmoModelGraph.ModelNode)
    indexdict = node.basenode.indices
    length(indexdict) > 1 && error("More than one index found for node")
    return collect(values(node.basenode.indices))[1]
end

function getgraph(node::Plasmo.PlasmoModelGraph.ModelNode)
    indexdict = node.basenode.indices
    length(indexdict) > 1 && error("More than one index found for node")
    return collect(keys(node.basenode.indices))[1]
end

function subgraphobjective(node, graph)
    if out_degree(graph, node) == 0
        return getobjectivevalue(getmodel(node))
    end
    return getobjectivevalue(getmodel(node)) + sum(subgraphobjective(childnode, graph) for childnode in out_neighbors(graph, node))
end
=#
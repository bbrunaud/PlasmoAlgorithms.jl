using PyCall
using Plasmo
using JuMP
using Gurobi

@pyimport smps.read as smps
pwd()

smodel = smps.StochasticModel("test/smps/sizes10")
numstages = length(smodel[:periods])
stages = 1:numstages
numvars = [length(smodel[:stage_vars][i]) for i in stages]
numconstrs = [length(smodel[:stage_constrs][i]) for i in stages]

g = PlasmoGraph()
g.attributes[:nodelabels] = String[]
ndl = Dict()

scenarios = smodel[:scenarios]
for scen in values(scenarios)
    genealogy = scen[:genealogy]
    for k in length(genealogy):-1:1
        sclabel = genealogy[k]
        if !in(sclabel,g.attributes[:nodelabels])
            ns = add_node!(g)
            ns.label = sclabel
            push!(g.attributes[:nodelabels],ns.label)
            ndl[sclabel] = ns

            if sclabel == "ROOT"
                sc = smodel
                m = _scenariomodel(sc,root=true)
            else
                sc = scenarios[sclabel]
                m = _scenariomodel(sc)
            end
            m.solver = GurobiSolver()
            setmodel(ns,m)
        end

        if k < length(genealogy)
            n1 = ndl[genealogy[k]]
            n2 = ndl[genealogy[k+1]]
            add_edge(g,n1,n2)
            @linkconstraint(g, [i in keys(n1[:y])], n1[:y][i] == n2[:x][i])
        end
    end
end

g.solver = GurobiSolver()

function _getstage(scenario,scmodel,root=false)
    root && return 1
    stage = 0
    stagename = scenario[:branch_period]
    stagenames = scmodel[:periods]
    stage = find(x -> x .== stagename, stagenames)[1]
end

function _getmatrix(stagedict::Dict,I,J)
    A = sparse(zeros(I,J))
    for (i,j) in keys(stagedict)
        A[i+1,j+1] = stagedict[i,j][1]
    end
    A
end

function _scenariomodel(sc;root=false)
    m = Model()
    stage = _getstage(sc,smodel,root)
    @variable(m, y[1:numvars[stage]])
    setlowerbound.(y,sc[:lb][stage])
    setupperbound.(y,sc[:ub][stage])
    for i in 1:numvars[stage]
        category = smodel[:vtype][stage][i]
        if category == "B"
            setcategory(y[i],:Bin)
        elseif category == "I"
            setcategory(y[i],:Int)
        end
    end
    c2 = sc[:c][stage]
    probability = root ? 1 : sc[:probability]
    @objective(m, Min, probability*c2'*y)

    A2 = _getmatrix(sc[:A][stage,stage],numconstrs[stage],numvars[stage])
    b = sc[:b][stage]

    if stage > 1
        @variable(m, x[1:numvars[stage-1]])
        setlowerbound.(x,sc[:lb][stage-1])
        setupperbound.(x,sc[:ub][stage-1])
        for i in 1:numvars[stage-1]
            category = smodel[:vtype][stage-1][i]
            if category == "B"
                setcategory(x[i],:Bin)
            elseif category == "I"
                setcategory(x[i],:Int)
            end
        end

        A1 = _getmatrix(sc[:A][stage,stage-1],numconstrs[stage],numvars[stage-1])
        @constraint(m, A1*x + A2*y .<= b)
    else
        @constraint(m, A2*y .<= b)
    end

    for k in 1:numconstrs[stage]
        c = m.linconstr[k]
        sense = sc[:b_sense][stage][k][1]
        if sense == ">"
            c.lb = c.ub
            c.ub = +Inf
        elseif sense == "="
            c.lb = c.ub
        end
    end

    return m
end

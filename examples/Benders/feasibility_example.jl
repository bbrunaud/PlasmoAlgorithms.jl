using JuMP
using Gurobi
using Plasmo
using PlasmoAlgorithms

graph = OptiGraph()
optimizer = Gurobi.Optimizer
const env = Gurobi.Env()
Plasmo.set_optimizer(graph, with_optimizer(optimizer, env))

mp = @optinode(graph)
JuMP.set_optimizer(mp.model, with_optimizer(optimizer, env))
set_optimizer_attribute(mp.model, "InfUnbdInfo", 1)
@variable(mp, y <= 5)
@objective(mp, Max, y)

sp = @optinode(graph)
JuMP.set_optimizer(sp.model, with_optimizer(optimizer, env))
set_optimizer_attribute(sp.model, "InfUnbdInfo", 1)
@variables(sp, begin
    v >= 1
    y >= 0
end)
@constraints(sp, begin
    v <= 3 - y
    v <= y + 3
end)
@objective(sp, Max, v)

edge1 = add_edge!(graph, [mp,sp])

@linkconstraint(graph, mp[:y] == sp[:y])

r = bendersoptimize!(graph, max_iterations = 5)

Gurobi.free_env(env)
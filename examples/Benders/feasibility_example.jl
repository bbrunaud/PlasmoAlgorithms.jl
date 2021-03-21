using JuMP
using Gurobi
using Plasmo
using PlasmoAlgorithms

graph = OptiGraph()
optimizer = Gurobi.Optimizer
Plasmo.set_optimizer(graph, optimizer)

mp = @optinode(graph)
JuMP.set_optimizer(mp.model, optimizer)
set_optimizer_attribute(mp, "InfUndbRay", 1)
@variable(mp, y <= 5)
@objective(mp, Max, y)

sp = @optinode(graph)
JuMP.set_optimizer(sp.model, optimizer)
set_optimizer_attribute(sp, "InfUndbRay", 1)
@variables(sp, begin
    v >= 1
    y >= 0
end)
@constraints(sp, begin
    v <= 3 - y
    v <= y + 3
end)
@objective(sp, Max, v)

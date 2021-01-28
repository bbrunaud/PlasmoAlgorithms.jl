using JuMP
using Xpress
using Plasmo
using PlasmoAlgorithms

graph = OptiGraph()
optimizer = Xpress.Optimizer
Plasmo.set_optimizer(graph, optimizer)


mp = @optinode(graph)
JuMP.set_optimizer(mp.model, optimizer)
@variable(mp,x[1:2]>=0)
@constraint(mp,x[1]+x[2]>=27)
@objective(mp,Min, 3x[1]+2x[2])

sp1 = @optinode(graph)
JuMP.set_optimizer(sp1.model, optimizer)
@variable(sp1,x[1:2]>=0)
@variable(sp1,y[1:2]>=0)
@constraint(sp1,x[1]+y[2]>=34)
@constraint(sp1,x[2]+y[1]>=52)
@objective(sp1,Min, 4y[1]+6y[2])

sp2 = @optinode(graph)
JuMP.set_optimizer(sp2.model, optimizer)
@variable(sp2,x[1:2]>=0)
@variable(sp2,z[1:2]>=0)
@constraint(sp2,x[1]+z[2]>=61)
@constraint(sp2,x[2]+z[1]>=29)
@objective(sp2,Min, 7z[1]+8z[2])

sp3 = @optinode(graph)
JuMP.set_optimizer(sp3.model, optimizer)
@variable(sp3,x[1:2]>=0)
@variable(sp3,t[1:2]>=0)
@constraint(sp3,x[1]+t[2]>=43)
@constraint(sp3,x[2]+t[1]>=75)
@objective(sp3,Min, 10t[1]+5t[2])

edge1 = add_edge!(graph, [mp,sp1])
edge2 = add_edge!(graph, [mp,sp2])
edge3 = add_edge!(graph, [mp,sp3])

@linkconstraint(graph, [i in 1:2], mp[:x][i] == sp1[:x][i])
@linkconstraint(graph, [i in 1:2], mp[:x][i] == sp2[:x][i])
@linkconstraint(graph, [i in 1:2], mp[:x][i] == sp3[:x][i])

r = bendersoptimize!(graph, max_iterations = 10)
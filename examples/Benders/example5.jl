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
@constraint(mp, x[1] + x[2] >= 50)
@objective(mp, Min, x[1] + 2x[2])

sp1 = @optinode(graph)
JuMP.set_optimizer(sp1.model, optimizer)
@variable(sp1,x[1:2]>=0)
@variable(sp1,y[1:2]>=0)
@constraint(sp1, x[1]+y[1]>=60)
@constraint(sp1, x[2]+y[2]>=70)
@objective(sp1, Min, 3y[1]+4y[2])

sp2 = @optinode(graph)
JuMP.set_optimizer(sp2.model, optimizer)
@variable(sp2,y[1:2]>=0)
@variable(sp2, z[1:2]>=0)
@constraint(sp2, y[1]+z[1]>=80)
@constraint(sp2, y[2]+z[2]>=90)
@objective(sp2, Min, 5z[1]+6z[2])

sp3 = @optinode(graph)
JuMP.set_optimizer(sp3.model, optimizer)
@variable(sp3,x[1:2]>=0)
@variable(sp3,t[1:2]>=0)
@constraint(sp3,x[1]+t[2]>=43)
@constraint(sp3,x[2]+t[1]>=75)
@objective(sp3,Min,.5t[1]+5t[2])


edge1 = add_edge!(graph, [mp,sp1])
edge2 = add_edge!(graph, [mp,sp3])
edge3 = add_edge!(graph, [sp1,sp2])

@linkconstraint(graph, [i in 1:2], mp[:x][i] == sp1[:x][i])
@linkconstraint(graph, [i in 1:2], mp[:x][i] == sp3[:x][i])
@linkconstraint(graph, [i in 1:2], sp1[:y][i] == sp2[:y][i])

r = bendersoptimize!(graph, max_iterations = 10)
using JuMP
using Xpress
using Plasmo
using PlasmoAlgorithms

graph = OptiGraph()
optimizer = Xpress.Optimizer
Plasmo.set_optimizer(graph, optimizer)


mp = @optinode(graph)
JuMP.set_optimizer(mp.model, optimizer)
@variable(mp,x[1:3]>=0)
@constraint(mp,x[1]+x[2]+x[3]<=500)
@objective(mp,Min,150x[1]+230x[2]+260x[3])

sp1 = @optinode(graph)
JuMP.set_optimizer(sp1.model, optimizer)
@variable(sp1, w[1:4]>=0)
@variable(sp1, y[1:2]>=0)
@variable(sp1, x[1:3]>=0)
@constraint(sp1, 3x[1]+y[1]-w[1]>=200)
@constraint(sp1, 3.6x[2]+y[2]-w[2]>= 240)
@constraint(sp1, w[3]+w[4]<=24x[3])
@constraint(sp1, w[3]<=6000)
@objective(sp1,Min,(-1/3)*(170w[1]-238y[1]+150w[2]-210y[2]+36w[3]+10w[4]))

sp2 = @optinode(graph)
JuMP.set_optimizer(sp2.model, optimizer)
@variable(sp2, w[5:8]>=0)
@variable(sp2, y[3:4]>=0)
@variable(sp2, x[1:3]>=0)
@constraint(sp2, 2.5x[1]+y[3]-w[5]>=200)
@constraint(sp2, 3x[2]+y[4]-w[6]>=240)
@constraint(sp2, w[7]+w[8]<=20x[3])
@constraint(sp2, w[7]<=6000)
@objective(sp2, Min, (-1/3)*(170w[5]-238y[3]+150w[6]-210y[4]+36w[7]+10w[8]))

sp3 = @optinode(graph)
JuMP.set_optimizer(sp3.model, optimizer)
@variable(sp3,w[9:12]>=0)
@variable(sp3,y[5:6]>=0)
@variable(sp3, x[1:3]>=0)
@constraint(sp3, 2x[1]+y[5]-w[9]>=200)
@constraint(sp3, 2.4x[2]+y[6]-w[10]>=240)
@constraint(sp3, w[11]+w[12]<=16x[3])
@constraint(sp3, w[11]<=6000)
@objective(sp3, Min, (-1/3)*(170w[9]-238y[5]+150w[10]-210y[6]+36w[11]+10w[12]))

edge1 = add_edge!(graph, [mp,sp1])
edge2 = add_edge!(graph, [mp,sp2])
edge3 = add_edge!(graph, [mp,sp3])


@linkconstraint(graph, [i in 1:3], mp[:x][i] == sp1[:x][i])
@linkconstraint(graph, [i in 1:3], mp[:x][i] == sp2[:x][i])
@linkconstraint(graph, [i in 1:3], mp[:x][i] == sp3[:x][i])

r = bendersoptimize!(graph, max_iterations = 10)
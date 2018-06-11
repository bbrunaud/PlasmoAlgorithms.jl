#Example 8

using JuMP
using Gurobi
using Plasmo
using PlasmoAlgorithms

mp = Model(solver = GurobiSolver())
sp1 = Model(solver = GurobiSolver())
sp2 = Model(solver = GurobiSolver())
sp3 = Model(solver = GurobiSolver())

@variable(mp, x >=0)
@objective(mp, Min, x)

@variable(sp1, y>=0)
@variable(sp1, x>=0)
@constraint(sp1, x + y >= 15)
@objective(sp1, Min, 2y)

@variable(sp2, x>=0)
@variable(sp2, z>=0)
@constraint(sp2, x + z <= 13)
@objective(sp2, Min, 3z)

@variable(sp3, z>=0)
@variable(sp3, t>=0)
@constraint(sp3, z + t<=3)
@objective(sp3, Min, 4t)

g = PlasmoGraph()
g.solver = GurobiSolver()
n1 = add_node(g)
n2 = add_node(g)
n3 = add_node(g)
n4 = add_node(g)

setmodel(n1, mp)
setmodel(n2, sp1)
setmodel(n3, sp2)
setmodel(n4, sp3)

@linkconstraint(g, n1[:x] == n2[:x])
@linkconstraint(g, n1[:x] == n3[:x])
@linkconstraint(g, n3[:z] == n4[:z])

edge = Plasmo.add_edge(g,n1,n2)
edge = Plasmo.add_edge(g,n1,n3)
edge = Plasmo.add_edge(g,n3,n4)

bendersolve(g)

using JuMP
using Gurobi
using Plasmo
using PlasmoAlgorithms

mp = Model(solver = GurobiSolver())
sp1 = Model(solver = GurobiSolver())
sp2 = Model(solver = GurobiSolver())

@variable(mp, x[1:2]>=0)
@objective(mp, Min, x[1]+2x[2])

@variable(sp1, x[1:2]>=0)
@variable(sp1, y>=0)
@constraint(sp1, x[1]+x[2]+y>=6)
@objective(sp1, Min, 3y)

@variable(sp2, x[1:2]>=0)
@variable(sp2, w>=0)
@constraint(sp2, -3x[1]+2x[2]+w>=7)
@objective(sp2, Min, 4w)

g = PlasmoGraph()
g.solver = GurobiSolver()

n1 = add_node(g)
n2 = add_node(g)
n3 = add_node(g)

setmodel(n1, mp)
setmodel(n2, sp1)
setmodel(n3, sp2)

edge = Plasmo.add_edge(g, n1, n2)
edge = Plasmo.add_edge(g, n1, n3)

@linkconstraint(g, [i in 1:2], n1[:x][i] == n2[:x][i])
@linkconstraint(g, [i in 1:2], n1[:x][i] == n3[:x][i])

test = lagrangePrep(g, [0,0])
#test = crossSolve(g,5)

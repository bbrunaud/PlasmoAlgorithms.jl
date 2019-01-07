using JuMP
using Gurobi
using Plasmo
using PlasmoAlgorithms

##Place MP and SP into ModelGraph
mp = Model(solver = GurobiSolver())
sp = Model(solver = GurobiSolver())

@variable(mp, x[1:2], category = :Bin)
@constraint(mp, x[1]+x[2]>=1)
@objective(mp, Min, 2x[1]+3x[2])

@variable(sp, y[1:2], category = :Bin)
@variable(sp, x[1:2], category = :Bin)
@constraint(sp, x[1]+x[2]+y[1]>=2)
@constraint(sp, y[1]+y[2]<=2)
@objective(sp, Min, y[1]+y[2])

## Plasmo Graph
g = ModelGraph()
setsolver(g, GurobiSolver())
n1 = add_node(g)
setmodel(n1,mp)
n2 = add_node(g)
setmodel(n2,sp)

##Set n2 as a child node of n1
edge = Plasmo.add_edge(g,n1,n2)

## Linking constraints between MP and SP
@linkconstraint(g,[i in 1:2], n1[:x][i] == n2[:x][i])

r = bendersolve(g)

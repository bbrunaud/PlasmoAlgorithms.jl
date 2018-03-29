using JuMP
using Gurobi
using Plasmo
using PlasmoAlgorithms
using Logging

Logging.configure(level=DEBUG)

##Place MP and SP into PlasmoGraph
mp = Model(solver = GurobiSolver())
sp = Model(solver = GurobiSolver())

@variable(mp,y>=0)
@objective(mp,Min,2y)

@variable(sp,x[1:2]>=0)
@variable(sp,y>=0)
@variable(sp, z[1:2], Bin)
@constraint(sp,c1,2x[1]-x[2]+3y>=4)
@constraint(sp,x[1]+2x[2]+y>=3)
@constraint(sp, bm[i=1:2], x[i] <= 10z[i])
@constraint(sp, one[i=1:2], z[1] + z[2] <= 1)
@objective(sp,Min,2x[1]+3x[2])


## Plasmo Graph
g = PlasmoGraph()
g.solver = GurobiSolver()
n1 = add_node(g)
setmodel(n1,mp)
n2 = add_node(g)
setmodel(n2,sp)

##Set n2 as a child node of n1
edge = Plasmo.add_edge(g,n1,n2)

## Linking constraints between MP and SP
@linkconstraint(g, n1[:y] == n2[:y])

bendersolve(g, max_iterations=5)


nlp = Model(solver = IpoptSolver())

@variable(nlp,x[1:2]>=0, upperbound=10)
@variable(nlp,y>=0)
@variable(nlp, z[1:2] >=0, upperbound=1)
@variable(nlp, yt >=0)
@constraint(nlp,c1,2x[1]-x[2]+3y>=4)
@constraint(nlp,x[1]+2x[2]+y>=3)
@NLconstraint(nlp, bm[i=1:2], 1-exp(-x[i]/0.01) <= z[i])
@constraint(nlp, one[i=1:2], z[1] + z[2] <= 1)
@constraint(nlp, dualcon, yt - y == 0)
@objective(nlp,Min,2x[1]+3x[2])

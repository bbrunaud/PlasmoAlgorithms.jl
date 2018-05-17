using JuMP
using Gurobi
using Plasmo

## Model from Fisher,1985. An Applications Oriented Guide to Lagrangian Relaxation
# Max 16x[1] + 10x[2] + 4y[2]
# s.t. x[1] + x[2] <= 1
#      y[1] + y[2] <= 1
#      8x[1] + 2x[2] + y[2] + 4y[2] <= 10
#      x, y ∈ {0,1}

## Model on x
# Min 16x[1] + 10x[2]
# s.t. x[1] + x[2] <= 1OutputFlag=0
#      x ∈ {0,1}
m1 = Model(solver=GurobiSolver(OutputFlag=0))

@variable(m1, xm[i in 1:2],Bin)
@constraint(m1, xm[1] + xm[2] <= 1)
@objective(m1, Max, 16xm[1] + 10xm[2])

## Model on y`
# Max  4y[2]
# s.t. y[1] + y[2] <= 1
#      8x[1] + 2x[2] + y[2] + 4y[2] <= 10
#      x, y ∈ {0,1}

m = Model(solver=GurobiSolver(OutputFlag=0))

@variable(m, xs[i in 1:2],Bin)
@variable(m, y[i in 1:2], Bin)
@constraint(m, y[1] + y[2] <= 1)
@constraint(m, 8xs[1] + 2xs[2] + y[2] + 4y[2] <= 10)
@objective(m, Max, 4y[2])

## Plasmo Graph
g = PlasmoGraph()
g.solver = GurobiSolver(OutputFlag=0)
n1 = add_node(g)
setmodel(n1,m1)
n2 = add_node(g)
setmodel(n2,m)


## Linking
# m1[x] = m2[x]  ∀i ∈ {1,2}
@linkconstraint(g, [i in 1:2], n1[:xm][i] == n2[:xs][i])

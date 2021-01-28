using Plasmo
using Xpress
using JuMP

graph = OptiGraph()
optimizer = Xpress.Optimizer

## Model on x
# Min 16x[1] + 10x[2]
# s.t. x[1] + x[2] <= 1OutputFlag=0
#      x ∈ {0,1}
n1 = @optinode(graph)
JuMP.set_optimizer(getmodel(n1), optimizer)
#set_silent(getmodel(n1))
@variable(n1, xm[i in 1:2],Bin)
@constraint(n1, xm[1] + xm[2] <= 1)
@objective(n1, Max, 16xm[1] + 10xm[2])

## Model on y`
# Max  4y[2]
# s.t. y[1] + y[2] <= 1
#      8x[1] + 2x[2] + y[2] + 4y[2] <= 10
#      x, y ∈ {0,1}
n2 = @optinode(graph)
JuMP.set_optimizer(getmodel(n2), optimizer)
#set_silent(getmodel(n2))
@variable(n2, xs[i in 1:2],Bin)
@variable(n2, y[i in 1:2], Bin)
@constraint(n2, y[1] + y[2] <= 1)
@constraint(n2, 8xs[1] + 2xs[2] + y[1] + 4y[2] <= 10)
@objective(n2, Max, 4y[2])

## Linking
# m1[x] = m2[x]  ∀i ∈ {1,2}
@linkconstraint(graph, [i in 1:2], n1[:xm][i] == n2[:xs][i])
Plasmo.set_optimizer(graph, optimizer)

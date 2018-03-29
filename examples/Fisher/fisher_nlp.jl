using JuMP
using BARON

m = Model(solver=BaronSolver())

@variable(m, xs[i in 1:2],Bin)
@variable(m, xm[i in 1:2],Bin)
@variable(m, λ[i in 1:2] >= -8, upperbound=8)
@variable(m, y[i in 1:2], Bin)
@variable(m, η)
@constraint(m, y[1] + y[2] <= 1)
@constraint(m, 8xs[1] + 2xs[2] + y[2] + 4y[2] <= 10)
@NLconstraint(m, η >= 4y[2] + 16xm[1] + 10xm[2] + sum(λ[i]*(xm[i] -xs[i]) for i in 1:2))
@objective(m, Min, η)

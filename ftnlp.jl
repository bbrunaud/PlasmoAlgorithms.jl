
include("ftbase.jl")

nlp = Model(solver=GurobiSolver())

@variable(nlp, x[1:3] >= 0)
@variable(nlp, y[1:3], >= 0, upperbound=1)
@variable(nlp, B >= 0, upperbound=600)
@constraint(nlp, sum(F[i]*y[i] for i in 1:3) <= B)
@NLconstraint(nlp, bm[i = 1:3], 1 - exp(-x[i]/1) == y[i])
@constraint(nlp, dem[i = 1:3], x[i] <= D[i])
@constraint(nlp, atmost2, sum(y[i] for i in 1:3) <= 2)
@objective(nlp, Min, B + sum(F[i]*y[i] - P*x[i] for i in 1:3))

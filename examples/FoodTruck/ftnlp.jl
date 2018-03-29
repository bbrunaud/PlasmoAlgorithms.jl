using Ipopt

include("ftbase.jl")

nlp = Model(solver=IpoptSolver())
@variable(nlp, x[i = 1:3] >= 0, upperbound=D[i])
@variable(nlp, y[1:3] >= 0, upperbound=1)
@variable(nlp, nBM >= 0, upperbound=600)
@variable(nlp, nB >= 0, upperbound=600)
@constraint(nlp, sum(F[i]*y[i] for i in 1:3) <= nB)
@NLconstraint(nlp, bm[i = 1:3], 1 - exp(-x[i]/1) == y[i])
@constraint(nlp, dem[i = 1:3], x[i] <= D[i])
@constraint(nlp, atmost2, sum(y[i] for i in 1:3) <= 2)
@objective(nlp, Min, nB + sum(F[i]*y[i] - P*x[i] for i in 1:3))

@constraint(nlp, dual, nBM == nB)

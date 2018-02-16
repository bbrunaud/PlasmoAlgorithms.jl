using JuMP
using Gurobi

F = [500,300,200]
D = [200,160,120]
B = 600
M = 200
P = 10

m = Model(solver=GurobiSolver())

@variable(m, x[1:3] >= 0)
@variable(m, y[1:3], Bin)
@variable(m, B >= 0, upperbound=600)
@constraint(m, sum(F[i]*y[i] for i in 1:3) <= B)
@constraint(m, bm[i = 1:3], x[i] <= M*y[i])
@constraint(m, dem[i = 1:3], x[i] <= D[i])
@objective(m, Min, B + sum(F[i]*y[i] - P*x[i] for i in 1:3))

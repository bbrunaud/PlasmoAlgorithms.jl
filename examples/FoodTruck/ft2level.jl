using JuMP
using Gurobi
using Ipopt

F = [500,300,200]
D = [200,160,120]
M = 200
P = 10

function opmodel(bval, modtype=:MILP)
m = Model()
@variable(m, x[1:3] >= 0)
if modtype == :NLP
    @variable(m, y[1:3] >= 0, upperbound=1)
else
    @variable(m, y[1:3], Bin)
end

@variable(m, b)
@variable(m, B)
@constraint(m, sum(F[i]*y[i] for i in 1:3) <= B)
@constraint(m, dem[i = 1:3], x[i] <= D[i])
@constraint(m, dl, B - b == 0)
@objective(m, Min, sum(F[i]*y[i] - P*x[i] for i in 1:3))

if modtype == :MILP
    @constraint(m, bm[i = 1:3], x[i] <= M*y[i])
    m.solver = GurobiSolver()
elseif modtype == :NLP
    @NLconstraint(m, bm[i = 1:3], exp(-x[i]/100) >= 1 - y[i])
    m.solver = IpoptSolver()
else
    error("Wrong model type")
end

setupperbound(b, bval)
setlowerbound(b, bval)
return m
end

ms = Model(solver=GurobiSolver())
@variable(ms, b >= 0, upperbound = 600)
@variable(ms, θ >= -3800)
@objective(ms, Min, θ)

function cut(ms, λ, bval, θk)
    θ = getindex(ms, :θ)
    b = getindex(ms, :b)
    return @LinearConstraint(θ >= θk - λ*(bval - b))
end

dual(m) = getdual(getindex(m,:dl))

using JuMP, Gurobi, OhMyREPL

mas = Model(solver=GurobiSolver())
@variable(mas, B >= 0, upperbound=600)
@variable(mas, θ >= -10000)
@objective(mas, Min, θ + 1.5B)

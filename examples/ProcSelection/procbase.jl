using JuMP
using Gurobi

FC = [1500 1000 3500]
Dem = 1

m = Model(solver=GurobiSolver())

@variable(m, B >= 0)
@variable(m, B2 >= 0, upperbound=0.8)
@variable(m, B3 >= 0, upperbound=0.8)
@variable(m, BP >= 0)
@variable(m, A >= 0)
@variable(m, A2 >= 0)
@variable(m, A3 >= 0)
@variable(m, C >= 0)
@variable(m, z[1:3], Bin)


@constraint(m, A == A2 + A3)
@constraint(m, B == B2 + B3 + BP)
@constraint(m, C <= Dem)

@constraint(m, B2 == 0.50A2)
@constraint(m, B3 == 0.70A3)
@constraint(m, C == 0.9B)

@constraint(m, C <= Dem*z[1])
@constraint(m, B2 <= 1*z[2])
@constraint(m, B3 <= 1*z[3])

@expression(m, θ,  - 13000C
                   + 3500z[1] + 2000C
                   + 1000z[2] + 1000B2
                   + 1500z[3] + 1200B3
                   + 1800A)
@objective(m, Min, 7000BP + θ)

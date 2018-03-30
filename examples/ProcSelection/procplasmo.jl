using JuMP
using Gurobi
using Plasmo
using PlasmoAlgorithms


FC = [1500 1000 3500]
Dem = 1

m1 = Model(solver=GurobiSolver())
@variable(m1, BPurch >= 0)
@objective(m1, Min, 7000BPurch)

m2 = Model(solver=GurobiSolver())

@variable(m2, B >= 0)
@variable(m2, B2 >= 0, upperbound=0.8)
@variable(m2, B3 >= 0, upperbound=0.8)
@variable(m2, BP >= 0)
@variable(m2, A >= 0)
@variable(m2, A2 >= 0)
@variable(m2, A3 >= 0)
@variable(m2, C >= 0)
@variable(m2, z[1:3], Bin)


@constraint(m2, A == A2 + A3)
@constraint(m2, B == B2 + B3 + BP)
@constraint(m2, C <= Dem)

@constraint(m2, B2 == 0.50A2)
@constraint(m2, B3 == 0.70A3)
@constraint(m2, C == 0.9B)

@constraint(m2, C <= Dem*z[1])
@constraint(m2, B2 <= 1*z[2])
@constraint(m2, B3 <= 1*z[3])

@objective(m2, Min, - 13000C
                   + 3500z[1] + 2000C
                   + 1000z[2] + 1000B2
                   + 1500z[3] + 1200B3
                   + 1800A)

g = PlasmoGraph()
n1 = add_node(g,m1)
n2 = add_node(g,m2)
add_edge(g,n1,n2)

@linkconstraint(g, n1[:BPurch] == n2[:BP])

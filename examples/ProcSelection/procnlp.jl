using Ipopt

m3 = Model(solver=IpoptSolver())

@variable(m3, B >= 0)
@variable(m3, B2 >= 0, upperbound=0.8)
@variable(m3, B3 >= 0, upperbound=0.8)
@variable(m3, BP >= 0)
@variable(m3, A >= 0)
@variable(m3, A2 >= 0)
@variable(m3, A3 >= 0)
@variable(m3, C >= 0, upperbound=Dem)
@variable(m3, z[1:3]>=0, upperbound=1)
@variable(m3, valbar >= 0)


@constraint(m3, A == A2 + A3)
@constraint(m3, B == B2 + B3 + BP)
@constraint(m3, C <= Dem)

@constraint(m3, B2 == 0.50A2)
@constraint(m3, B3 == 0.70A3)
@constraint(m3, C == 0.9B)

@NLconstraint(m3, 1 - exp(-C/0.1) <= z[1])
@NLconstraint(m3, 1 - exp(-B2/0.1) <= z[2])
@NLconstraint(m3, 1 - exp(-B3/0.1) <= z[3])

@constraint(m3, dualcon, valbar - BP == 0)

@objective(m3, Min, - 13000C
                   + 3500z[1] + 2000C
                   + 1000z[2] + 1000B2
                   + 1500z[3] + 1200B3
                   + 1800A)

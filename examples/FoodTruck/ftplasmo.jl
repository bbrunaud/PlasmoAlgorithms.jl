using Plasmo
using PlasmoAlgorithms

include("ftbase.jl")

m1 = Model(solver=GurobiSolver())
@variable(m1, BM >= 0, upperbound=600)
@objective(m1, Min, 1.5BM)

@objective(m, Min, sum(F[i]*y[i] - P*x[i] for i in 1:3))

m2 = deepcopy(m)

g = PlasmoGraph()
n1 = add_node(g, m1)
n2 = add_node(g, m)

add_edge(g, n1,n2)

@linkconstraint(g, n1[:BM] == n2[:B])

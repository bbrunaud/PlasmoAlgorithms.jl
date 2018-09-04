using PlasmoAlgorithms
using Base.Test

const PA = PlasmoAlgorithms

include("fisher_min.jl")

c = PA.CrossGraph(g)
@test typeof(c) == PA.CrossGraph
PA.crossprepare(c)

bn1 = getnode(c.bd,1)
cn1 = getnode(c.lg, 1)

# Begin Iterations
for i in 1:5
    ls = lagrangesolve(c.lg, max_iterations=1, initialmultipliers=:relaxation, lagrangeheuristic=x->-16)
    print(ls.termination)
    PA.lagrange_to_benders(c, termination=ls.termination)
end

bs = bendersolve(c.bd, max_iterations=1)
@test bs.bestbound == -16

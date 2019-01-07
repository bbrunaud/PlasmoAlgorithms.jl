using PlasmoAlgorithms
using Base.Test

const PA = PlasmoAlgorithms

include("example5.jl")

c = PA.CrossGraph(g)
@test typeof(c) == PA.CrossGraph
PA.crossprepare(c)

bn1 = getnode(c.bd,1)
cn1 = getnode(c.lg, 1)

# Begin Iterations
PA.bendersolve(c.bd, max_iterations=1)
@test length(getattribute(bn1, :cutdata)[2]) == 1
PA.benders_to_lagrange(c)
@test length(getattribute(c.lg, :cutdata)[2]) == 1
PA.lagrangesolve(c.lg, max_iterations=1)
@test length(getattribute(c.lg, :cutdata)[2]) == 2
PA.lagrange_to_benders(c)

# 2 cuts from Benders and 2 cuts from Lagrange
@test length(getattribute(c.lg, :cutdata)) == 4


using PlasmoAlgorithms
using Base.Test

include("benders_example.jl")

r = bendersolve(g)

@test isapprox(r.objval, 5.28, atol=0.01)

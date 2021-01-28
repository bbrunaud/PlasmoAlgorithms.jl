using PlasmoAlgorithms
using Test

include("fisher.jl")

δ = 0.5
maxnoimprove = 3
cpbound=1e6

heur(g) = 16
sg = lagrangeoptimize!(graph, update_method=:subgradient, lagrangeheuristic=heur, verbose=false)

@test sg.objval == 16
@test sg.bestbound == 16
@test sg.gap == 0

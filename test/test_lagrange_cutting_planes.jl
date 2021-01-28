using PlasmoAlgorithms
using Test

include("fisher.jl")

δ = 0.5
maxnoimprove = 3
cpbound=1e6

include("fisher.jl")
cp = lagrangeoptimize!(graph, update_method=:cuttingplanes, lagrangeheuristic=heur, verbose=false)

@test sg.objval == 16
@test sg.bestbound == 16
@test sg.gap == 0

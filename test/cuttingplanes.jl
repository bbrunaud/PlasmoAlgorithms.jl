using PlasmoAlgorithms

include("fisher.jl")

heur(g) = 16

r = lagrangesolve(g, update_method=:cuttingplanes, lagrangeheuristic=heur, cpbound=16)

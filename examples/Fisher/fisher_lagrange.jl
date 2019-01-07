using Plasmo
using JuMP
using PlasmoAlgorithms

include("fisher.jl")

heur(g) = 16

r = lagrangesolve(g,update_method=:subgradient,max_iterations=30,lagrangeheuristic=heur)

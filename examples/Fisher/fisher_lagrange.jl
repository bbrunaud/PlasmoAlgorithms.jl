@everywhere using Plasmo
@everywhere using PlasmoAlgorithms
@everywhere using JuMP
@everywhere using Logging

using OhMyREPL

include("fisher.jl")


Logging.configure(level=DEBUG)

heur(g) = 16
#r = lagrangesolve(g,update_method=:subgradient,max_iterations=10,lagrangeheuristic=heur)

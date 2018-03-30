using Plasmo
using JuMP
using Logging
using PlasmoAlgorithms

using OhMyREPL

include("fisher.jl")


Logging.configure(level=DEBUG)


heur(g) = 16
#r = lagrangesolve(g,update_method=:subgradient,max_iterations=30,lagrangeheuristic=heur)

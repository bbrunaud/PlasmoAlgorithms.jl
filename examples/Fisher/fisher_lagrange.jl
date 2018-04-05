@everywhere using Plasmo
@everywhere using JuMP
@everywhere using Logging

include("fisher.jl")


Logging.configure(level=DEBUG)

#r = lagrangesolve(g,update_method=:subgradient,max_iterations=10)

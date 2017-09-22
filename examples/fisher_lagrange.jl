@everywhere using Plasmo
@everywhere using JuMP
@everywhere using Logging

include("fisher.jl")


Logging.configure(level=DEBUG)

r = lagrangesolve(graph,update_method=:subgradient,max_iterations=10)
println(r)

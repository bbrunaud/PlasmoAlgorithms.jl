include("toy.jl")
include("utils.jl")
include("solution.jl")
include("benders.jl")



bendersolve(g, max_iterations=10, cuts=[:LP])
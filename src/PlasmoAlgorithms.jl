module PlasmoAlgorithms

using Plasmo
using JuMP
using Requires


import Plasmo.solve, Base.==

export Solution, lagrangesolve, bendersolve

include("types.jl")
include("lagrange.jl")
include("benders.jl")
include("solution.jl")
include("utils.jl")

end

module PlasmoAlgorithms


using Plasmo
using JuMP
using Requires
using PyCall
using CPLEX
using MathProgBase

import Plasmo.solve, Base.==

export Solution, lagrangesolve, psolve, bendersolve,
lgprepare, solvenode,crosssolve,fixbinary, nearest_scenario,random_scenario, bab_solve,BABnode,
#largest_rel_diameter,largest_distance,set_bounds,solve_node,copy_node,
calculate_gap,
# Solution
saveiteration,

# Utils
normalizegraph, smpsread,

#convexification
add_PiecewiseMcCormick, add_McCormick

include("lagrange.jl")
include("benders.jl")
include("cross.jl")
include("solution.jl")
include("utils.jl")
include("spatialbranchandbound.jl")
include("convexrelaxations.jl")

end

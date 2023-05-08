module PlasmoAlgorithms

    using Plasmo
    using JuMP
    using Requires
    using LinearAlgebra
    using Distributed
    using Plots

    import Base.==

    export Solution, lagrangeoptimize!, bendersoptimize!, crossoptimize!

    include("utils.jl")
    include("types.jl")
    include("lagrange.jl")
    include("benders.jl")
    include("cross.jl")
    include("solution.jl")


end

include("fisher_min.jl")

using PlasmoAlgorithms
const PA = PlasmoAlgorithms
const MPB = MathProgBase

PA.bdprepare(g)
θ = getindex(m1, :θ)
setlowerbound(θ[2],-40)
solve(m1)
LB,UB = PA.forwardstep(g)
PA.backwardstep(g,:Root)

LB,UB = PA.forwardstep(g)

PA.rootcut(g,n2)
LB,UB = PA.forwardstep(g)

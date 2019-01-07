using PlasmoAlgorithms
using Base.Test

const PA = PlasmoAlgorithms

include("fisher.jl")

δ = 0.5
maxnoimprove = 3
cpbound=1e6

@test typeof(g) == ModelGraph

node1 = getnode(g,1)
@test typeof(node1) == ModelNode

# Prepare Function
@test PA.lgprepare(g, δ, maxnoimprove, cpbound)
m1 = getmodel(node1)
@test m1.objSense == :Min
@test getattribute(g, :numlinks) == 2
@test typeof(m1.solver) == Gurobi.GurobiSolver
@test getattribute(g, :normalized) == -1

# Inital Relaxation
λ0 = getattribute(g,:λ)[end]
@test λ0 == zeros(2)
lpobj = PA.initialrelaxation(g)
@test lpobj == -18.0
λlp = getattribute(g,:λ)[end]
@test λlp == [8.0, 2.0]

nodes = [node for node in collect(getnodes(g))]
@test length(nodes) == 2

heur(g) = 16
r = lagrangesolve(g, update_method=:subgradient, lagrangeheuristic=heur)

@test r.objval == 16
@test r.bestbound == 16
@test r.gap == 0

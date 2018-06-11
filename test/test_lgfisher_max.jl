using PlasmoAlgorithms
using Base.Test

const PA = PlasmoAlgorithms

include("fisher.jl")

δ = 0.5
maxnoimprove = 3

node1 = getnodes(g)[1]

gc = deepcopy(g)

# Prepare Function
PA.lgprepare(g, δ, maxnoimprove)
@test node1.model.objSense == :Min
@test g.attributes[:numlinks] == 2
@test typeof(node1.model.solver) == Gurobi.GurobiSolver
@test g.attributes[:normalized] == -1

# Inital Relaxation
λ0 = g.attributes[:λ][end]
@test λ0 == zeros(2)
lpobj = PA.initialrelaxation(g)
@test lpobj == -17.6
λlp = g.attributes[:λ][end]
#@test λlp

nodes = [node for node in values(getnodes(g))]
@test length(nodes) == 2

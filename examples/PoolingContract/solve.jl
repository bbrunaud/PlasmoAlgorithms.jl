using PlasmoAlgorithms
using Plasmo
using CPLEX
using BARON
include("input.jl")
include("benderssub.jl")
include("bendersmaster.jl")
include("lagsub.jl")
include("ubsub.jl")


master = generate_bendersmaster()
g = PlasmoGraph()
g.solver = CplexSolver()
n1 = add_node(g)
setmodel(n1, master)

for s in scenarios
	benderssub = generate_benderssub(psi_f = psi_f, psi_d1=psi_d1, psi_d2=psi_d2, psi_b1=psi_b1, psi_b2=psi_b2, DU=DU[:, s], prob=prob[s])
	ubsub = generate_ubsub(psi_f = psi_f, psi_d1=psi_d1, psi_d2=psi_d2, psi_b1=psi_b1, psi_b2=psi_b2, DU=DU[:, s], prob=prob[s])
	n = add_node(g)
	setmodel(n, benderssub)
	add_edge(g, n1, n)
	A = getindex(n1.model,:A)
	@linkconstraint(g, [i in feeds], n1[:A][i] == n[:A][i])
	gamma_intlt = getindex(n1.model, :gamma_intlt)
	@linkconstraint(g, [i in feeds], n1[:gamma_intlt][i] == n[:gamma_intlt][i])
	gamma_pool = getindex(n1.model, :gamma_pool)
	@linkconstraint(g, [i in pools], n1[:gamma_pool][i] == n[:gamma_pool][i])
	S = getindex(n1.model, :S)
	@linkconstraint(g, [i in pools], n1[:S][i] == n[:S][i])
	n.attributes[:ubsub] = ubsub
end

bendersolve(g, cuts=[:LP], max_iterations=60, timelimit=100000, is_nonconvex=true)

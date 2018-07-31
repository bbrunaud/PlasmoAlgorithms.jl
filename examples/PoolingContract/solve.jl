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
	n.attributes[:scenario] = s
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

# bendersolve(g, cuts=[:LP], max_iterations=60, timelimit=100000, is_nonconvex=true)

lag_g = PlasmoGraph()
lag_g.solver = BaronSolver(maxtime=5e4, epsr= 1e-3, CplexLibName = "/opt/ibm/ILOG/CPLEX_Studio127/cplex/bin/x86-64_linux/libcplex1270.so")
lagsub1 = generate_lagsub(psi_f = psi_f, psi_d1=psi_d1, psi_d2=psi_d2, psi_b1=psi_b1, psi_b2=psi_b2, DU=DU[:, 1], prob=prob[1])
n1 = add_node(lag_g)
setmodel(n1, lagsub1)
n1.attributes[:prob] =prob[1]
n1.attributes[:scenario] = 1
for s in scenarios
	if s == length(scenarios)
		break
	end
	ss = s + 1
	n = add_node(lag_g)
	n.attributes[:prob] = prob[ss]
	n.attributes[:scenario] = ss 
	lagsub = generate_lagsub(psi_f = psi_f, psi_d1=psi_d1, psi_d2=psi_d2, psi_b1=psi_b1, psi_b2=psi_b2, DU=DU[:, ss], prob=prob[ss])
	setmodel(n, lagsub)
	@linkconstraint(lag_g, [i in feeds], n1[:A][i] == n[:A][i])
	@linkconstraint(lag_g, [i in feeds], n1[:gamma_intlt][i] == n[:gamma_intlt][i])
	@linkconstraint(lag_g, [i in pools], n1[:gamma_pool][i] == n[:gamma_pool][i])
	@linkconstraint(lag_g, [i in pools], n1[:S][i] == n[:S][i])
end
lagrangesolve(lag_g, max_iterations=2, lagrangeheuristic=:nearest_scenario,  maxnoimprove = 1)
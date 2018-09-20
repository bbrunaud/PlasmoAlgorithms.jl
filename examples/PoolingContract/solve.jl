using PlasmoAlgorithms
using Plasmo
using CPLEX
using BARON
using Gurobi
include("input27.jl")
include("benderssub.jl")
include("bendersmaster.jl")
include("lagsub.jl")
include("ubsub.jl")


master = generate_bendersmaster()
g = PlasmoGraph()
g.solver = CplexSolver(CPX_PARAM_SCRIND=0, CPXPARAM_Simplex_Tolerances_Feasibility=1e-9)
n1 = add_node(g)
setmodel(n1, master)
println(master)
println(master.colNames)
for s in scenarios
	benderssub = generate_benderssub(psi_f = psi_f[:,s], psi_d1=psi_d1[:,s], psi_d2=psi_d2[:,s], psi_b1=psi_b1[:,s], psi_b2=psi_b2[:,s], d = d[:,s], y_up = y_up[:, :, s], DU=DU[:, s], prob=prob[s])
	ubsub = generate_ubsub(psi_f = psi_f[:,s], psi_d1=psi_d1[:,s], psi_d2=psi_d2[:,s], psi_b1=psi_b1[:,s], psi_b2=psi_b2[:,s], d = d[:,s],  DU=DU[:, s], prob=prob[s])
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
	println(benderssub)
	println(ubsub)
end

# bendersolve(g, cuts=[:LIFT], max_iterations=100, rel_gap=1e-3, timelimit=100000, is_nonconvex=true)

lag_g = PlasmoGraph()
lag_g.solver = BaronSolver(maxtime=5e4, epsr= 1e-3, CplexLibName = "/opt/ibm/ILOG/CPLEX_Studio127/cplex/bin/x86-64_linux/libcplex1270.so")
lagsub1 = generate_lagsub(psi_f = psi_f[:,1], psi_d1=psi_d1[:,1], psi_d2=psi_d2[:,1], psi_b1=psi_b1[:,1], psi_b2=psi_b2[:,1], d = d[:,1], DU=DU[:, 1], prob=prob[1])
n1 = add_node(lag_g)
setmodel(n1, lagsub1)
n1.attributes[:prob] =prob[1]
n1.attributes[:scenario] = 1
println(lagsub1)
for s in scenarios
	if s == length(scenarios)
		break
	end
	ss = s + 1
	n = add_node(lag_g)
	n.attributes[:prob] = prob[ss]
	n.attributes[:scenario] = ss 
	lagsub = generate_lagsub(psi_f = psi_f[:,ss], psi_d1=psi_d1[:,ss], psi_d2=psi_d2[:,ss], psi_b1=psi_b1[:,ss], psi_b2=psi_b2[:,ss], d = d[:,ss], DU=DU[:, ss], prob=prob[ss])
	setmodel(n, lagsub)
	@linkconstraint(lag_g, [i in feeds], n1[:A][i] == n[:A][i])
	@linkconstraint(lag_g, [i in feeds], n1[:gamma_intlt][i] == n[:gamma_intlt][i])
	@linkconstraint(lag_g, [i in pools], n1[:gamma_pool][i] == n[:gamma_pool][i])
	@linkconstraint(lag_g, [i in pools], n1[:S][i] == n[:S][i])
	println(lagsub)
end
# lagrangesolve(lag_g, max_iterations=5, lagrangeheuristic=nearest_scenario,  maxnoimprove = 1)
# crosssolve(g, lag_g, max_iterations_lag=2, max_iterations_benders=1, is_nonconvex=true)

# all_α = [2 1.5 1 0.5]
# all_δ = [0.8 0.5 0.3]
# # all_α = [2 ]
# # all_δ = [0.8 ]
# gap = []
# for α in all_α
# 	for δ in all_δ
# 		gg = deepcopy(g)
# 		lagg = deepcopy(lag_g)
# 		solution = crosssolve(gg, lagg, max_iterations_lag=30, max_iterations_benders=60, is_nonconvex=true, lag_α = α, lag_δ = δ)
# 		push!(gap, solution.gap)
# 	end
# end
# solution = crosssolve(g, lag_g, max_iterations_lag=30, max_iterations_benders=60, is_nonconvex=true, lag_α = 1.5, lag_δ = 0.3)
# println(gap)

#start debug spatial bab
rootnode = BABnode(lag_g, g, -1e6, 1e6, [], [], [], [], :notfathomed)
start = time()
results = bab_solve(rootnode,  benders_cuts=[:LP], max_iterations_lag=10, max_iterations_benders=60, is_nonconvex=true, lag_α = 1.5, lag_δ = 0.3, rel_gap=1e-3)
# results = bab_solve(rootnode,  benders_cuts=[:LIFT], max_iterations_lag=10, max_iterations_benders=60, is_nonconvex=true, lag_α = 1.5, lag_δ = 0.3, rel_gap=1e-3)
# results = bab_solve(rootnode,  heuristic=:lagonly, benders_cuts=[:LP], max_iterations_lag=10, max_iterations_benders=60, is_nonconvex=true, lag_α = 1.5, lag_δ = 0.3, rel_gap=1e-3)
endd = time()
walltime = endd - start



















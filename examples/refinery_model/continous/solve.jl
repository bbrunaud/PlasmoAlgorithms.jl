using PlasmoAlgorithms
using Plasmo
include("input.jl")
include("lagsub.jl")
include("benderssub.jl")
include("master.jl")
include("ubsub.jl")

master = generate_master()
g = PlasmoGraph()
g.solver = CplexSolver(CPX_PARAM_SCRIND=0, CPXPARAM_Simplex_Tolerances_Feasibility=1e-9, CPX_PARAM_EPRHS=1e-9)
n1 = add_node(g)
setmodel(n1, master)
println(master)
for s in scenarios
	benderssub = generate_benderssub(prob=prob[s], Crude_yield_data = Crude_yield_data[:,:,s], Desulphurisation_cost=Desulphurisation_cost[:,s], Sulphur_2=Sulphur_2[:,s], Sulphur_GO_data= Sulphur_GO_data[:,s])
	ubsub = generate_ubsub(prob=prob[s], Crude_yield_data = Crude_yield_data[:,:,s], Desulphurisation_cost=Desulphurisation_cost[:,s], Sulphur_2=Sulphur_2[:,s], Sulphur_GO_data= Sulphur_GO_data[:,s])
	n = add_node(g)
	n.attributes[:scenario] = s
	setmodel(n, benderssub)
	add_edge(g, n1, n)
	@linkconstraint(g, [c in crudes], n1[:crudeQuantity][c] == n[:crudeQuantity][c])
	@linkconstraint(g, [c in crudes], n1[:pickCrude][c] == n[:pickCrude][c])
	n.attributes[:ubsub] = ubsub
	println(benderssub)
	println(ubsub)
end
# bendersolve(g, cuts=[:LIFT], max_iterations=100, rel_gap=1e-3, timelimit=100000, is_nonconvex=true)
lag_g = PlasmoGraph()
lag_g.solver = BaronSolver(maxtime=1e4, epsr= 1e-4, CplexLibName = "/opt/ibm/ILOG/CPLEX_Studio127/cplex/bin/x86-64_linux/libcplex1270.so")

lagsub1 = generate_lagsub(prob=prob[1], Crude_yield_data = Crude_yield_data[:,:,1], Desulphurisation_cost=Desulphurisation_cost[:,1], Sulphur_2=Sulphur_2[:,1], Sulphur_GO_data= Sulphur_GO_data[:,1])
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
	lagsub = generate_lagsub(prob=prob[ss], Crude_yield_data = Crude_yield_data[:,:,ss], Desulphurisation_cost=Desulphurisation_cost[:,ss], Sulphur_2=Sulphur_2[:,ss], Sulphur_GO_data= Sulphur_GO_data[:,ss])
	setmodel(n, lagsub)
	@linkconstraint(lag_g, [c in crudes], n1[:crudeQuantity][c] == n[:crudeQuantity][c])
	@linkconstraint(lag_g, [c in crudes], n1[:pickCrude][c] == n[:pickCrude][c])
	println(lagsub)
end

# lagrangesolve(lag_g, max_iterations=100, lagrangeheuristic=nearest_scenario,  maxnoimprove = 1, user_UB=-18000)



rootnode = BABnode(lag_g, g, -1e6, 1e6, [], [], [], [], :notfathomed)
start = time()
results = bab_solve(rootnode, time_limit=1e4, benders_cuts=[:LIFT], lagrangeheuristic=random_scenario, max_iterations_lag=20, max_iterations_benders=60, is_nonconvex=true, lag_α = 1.5, lag_δ = 0.3, rel_gap=1e-3, user_UB=-18000)
# results = bab_solve(rootnode, time_limit=1e4, benders_cuts=[:LP], lagrangeheuristic=random_scenario, max_iterations_lag=20, max_iterations_benders=60, is_nonconvex=true, lag_α = 1.5, lag_δ = 0.3, rel_gap=1e-3, user_UB=-18000)
# results = bab_solve(rootnode, heuristic=:lagonly, time_limit=1e4, lagrangeheuristic=random_scenario, benders_cuts=[:LP], max_iterations_lag=20, max_iterations_benders=1, is_nonconvex=true, lag_α = 1.5, lag_δ = 0.3, rel_gap=1e-3, user_UB=-18000)
endd = time()
walltime = endd - start



























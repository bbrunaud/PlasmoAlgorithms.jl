using PlasmoAlgorithms
using Plasmo
using CPLEX
using BARON
include("input.jl")
include("benderssub.jl")
include("master.jl")
include("lagsub.jl")
include("ubsub.jl")

master = generate_master()
g = PlasmoGraph()
g.solver = CplexSolver(CPX_PARAM_SCRIND=0, CPXPARAM_Simplex_Tolerances_Feasibility=1e-9)
n1 = add_node(g)
setmodel(n1, master)
println(master)
for s in scenarios
	benderssub = generate_benderssub(prob = prob[s], DemandPerDay = DemandPerDay[:, s], TotalDemandPerDay=TotalDemandPerDay[s], costPerTon_ub = costPerTon_ub[s], costPerTon_lb = costPerTon_lb[s], cycleTime_ub = cycleTime_ub[s], cycleTime_lb = cycleTime_lb[s] )
	ubsub = generate_ubsub(prob = prob[s], DemandPerDay = DemandPerDay[:, s], TotalDemandPerDay=TotalDemandPerDay[s], costPerTon_ub = costPerTon_ub[s], costPerTon_lb = costPerTon_lb[s], cycleTime_ub = cycleTime_ub[s], cycleTime_lb = cycleTime_lb[s] )
	n = add_node(g)
	n.attributes[:scenario] = s
	setmodel(n, benderssub)
	add_edge(g, n1, n)
	@linkconstraint(g, [p in products], n1[:productTankSize][p] == n[:productTankSize][p])	
	@linkconstraint(g, n1[:firstStageCost] == n[:firstStageCost])
	n.attributes[:ubsub] = ubsub
	println(benderssub)
	println(ubsub)
end
# bendersolve(g, cuts=[:LIFT], max_iterations=100, rel_gap=1e-3, timelimit=100000, is_nonconvex=true)

lag_g = PlasmoGraph()
lag_g.solver = BaronSolver(maxtime=5e4, epsr= 1e-3, CplexLibName = "/opt/ibm/ILOG/CPLEX_Studio127/cplex/bin/x86-64_linux/libcplex1270.so")
lagsub1 = generate_lagsub(prob = prob[1], DemandPerDay = DemandPerDay[:, 1], TotalDemandPerDay=TotalDemandPerDay[1], costPerTon_ub = costPerTon_ub[1], costPerTon_lb = costPerTon_lb[1], cycleTime_ub = cycleTime_ub[1], cycleTime_lb = cycleTime_lb[1] , all_prob = prob, all_TotalDemandPerDay = TotalDemandPerDay)
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
	lagsub = generate_lagsub(prob = prob[ss], DemandPerDay = DemandPerDay[:, ss], TotalDemandPerDay=TotalDemandPerDay[ss], costPerTon_ub = costPerTon_ub[ss], costPerTon_lb = costPerTon_lb[ss], cycleTime_ub = cycleTime_ub[ss], cycleTime_lb = cycleTime_lb[ss] , all_prob = prob, all_TotalDemandPerDay = TotalDemandPerDay)
	setmodel(n, lagsub)
	@linkconstraint(lag_g, [p in products], n1[:productTankSize][p] == n[:productTankSize][p])
	@linkconstraint(lag_g, n1[:firstStageCost] == n[:firstStageCost])
	println(lagsub)
end

# lagrangesolve(lag_g, max_iterations=30, lagrangeheuristic=nearest_scenario,  maxnoimprove = 1)
# all_α = [2 1.5 1 0.5]
# all_δ = [0.8 0.5 0.3]
# graphs = []
# gap = []
# # all_α = [2 ]
# # all_δ = [0.8 ]
# gap = []
# for α in all_α
# 	for δ in all_δ
# 		lag_gg = deepcopy(lag_g)
# 		solution = lagrangesolve(lag_gg, max_iterations=100, lagrangeheuristic=nearest_scenario,  maxnoimprove = 1)
# 		push!(graphs, lag_gg)
# 		push!(gap, solution.gap)
# 	end
# end
rootnode = BABnode(lag_g, g, -1e6, 1e6, [], [], [], [], :notfathomed)
start = time()
# results = bab_solve(rootnode, time_limit=1e4, benders_cuts=[:LIFT], max_iterations_lag=1, max_iterations_benders=1, is_nonconvex=true, lag_α = 1.5, lag_δ = 0.3, rel_gap=1e-3)
# results = bab_solve(rootnode, time_limit=1e4, benders_cuts=[:LP], max_iterations_lag=1, max_iterations_benders=1, is_nonconvex=true, lag_α = 1.5, lag_δ = 0.3, rel_gap=1e-3)
results = bab_solve(rootnode, heuristic=:lagonly, time_limit=1e4, benders_cuts=[:LP], max_iterations_lag=1, max_iterations_benders=1, is_nonconvex=true, lag_α = 1.5, lag_δ = 0.3, rel_gap=1e-3)
endd = time()
walltime = endd - start


























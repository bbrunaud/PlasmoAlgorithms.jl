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
	productTankSize = getindex(n1.model, :productTankSize)
	@linkconstraint(g, [p in products], n1[:productTankSize][p] == n[:productTankSize][p])
	firstStageCost = getindex(n1.model, :firstStageCost)
	@linkconstraint(g, n1[:firstStageCost] == n[:firstStageCost])
	
	n.attributes[:ubsub] = ubsub
	println(benderssub)
	println(ubsub)
end

model1 = create_flat_graph_model(g)
ProductionLength_upper = ones(3) * 3
master2 = generate_master()
g2= PlasmoGraph()
g2.solver = CplexSolver(CPX_PARAM_SCRIND=0, CPXPARAM_Simplex_Tolerances_Feasibility=1e-9)
n1 = add_node(g2)
setmodel(n1, master2)
println(master2)
for s in scenarios
	benderssub = generate_benderssub(prob = prob[s], DemandPerDay = DemandPerDay[:, s], TotalDemandPerDay=TotalDemandPerDay[s], costPerTon_ub = costPerTon_ub[s], costPerTon_lb = costPerTon_lb[s], cycleTime_ub = cycleTime_ub[s], cycleTime_lb = cycleTime_lb[s] )
	ubsub = generate_ubsub(prob = prob[s], DemandPerDay = DemandPerDay[:, s], TotalDemandPerDay=TotalDemandPerDay[s], costPerTon_ub = costPerTon_ub[s], costPerTon_lb = costPerTon_lb[s], cycleTime_ub = cycleTime_ub[s], cycleTime_lb = cycleTime_lb[s] )
	n = add_node(g2)
	n.attributes[:scenario] = s
	setmodel(n, benderssub)
	add_edge(g2, n1, n)
	productTankSize = getindex(n1.model, :productTankSize)
	@linkconstraint(g2, [p in products], n1[:productTankSize][p] == n[:productTankSize][p])
	firstStageCost = getindex(n1.model, :firstStageCost)
	@linkconstraint(g2, n1[:firstStageCost] == n[:firstStageCost])
	
	n.attributes[:ubsub] = ubsub
	println(benderssub)
	println(ubsub)
end
model2 = create_flat_graph_model(g2)
model1.solver = BaronSolver()
model2.solver= BaronSolver()
solve(model1)
solve(model2)
# bendersolve(g, cuts=[:LIFT], max_iterations=100, rel_gap=1e-3, timelimit=100000, is_nonconvex=true)






























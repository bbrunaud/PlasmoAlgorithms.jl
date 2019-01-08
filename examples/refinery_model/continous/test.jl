using Plasmo
using PlasmoAlgorithms

include("input.jl")
include("lagsub.jl")
include("benderssub.jl")
include("master.jl")
include("ubsub.jl")

master = generate_master()
g = PlasmoGraph()
g.solver = CplexSolver(CPX_PARAM_SCRIND=1, CPXPARAM_Simplex_Tolerances_Feasibility=1e-9, CPX_PARAM_EPRHS=1e-9)
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
	# println(ubsub)
end

fullspacemodel = create_flat_graph_model(g)
fullspacemodel.solver = g.solver
solve(fullspacemodel)
println(getobjectivevalue(fullspacemodel))





























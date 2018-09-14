using PlasmoAlgorithms
using Plasmo
include("input.jl")
include("lagsub.jl")

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

lagrangesolve(lag_g, max_iterations=100, lagrangeheuristic=nearest_scenario,  maxnoimprove = 1)































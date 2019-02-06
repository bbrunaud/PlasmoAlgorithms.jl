using JuMP
using CPLEX
using Plasmo
using PlasmoAlgorithms
include("input.jl")
include("master.jl")
include("subproblem.jl")
#define master problem in financeplanning 
master = generate_master()

g = ModelGraph()
setsolver(g, CplexSolver())
n1 = add_node(g)
setmodel(n1, master)

for s1 in base_scenarios
	n2 = add_node(g)
	sub = generate_sub(0, Return,  s1, false)
	setmodel(n2, sub)
	Plasmo.add_edge(g, n1, n2)
	@linkconstraint(g, [i in Investments], n1[:AmountInvestedCurrentStage][i] == n2[:AmountInvestedPreviousStage][i])
	for s2 in base_scenarios
		n3 = add_node(g)
		sub = generate_sub(0.25, Return,  s2, true)
		setmodel(n3, sub)
		Plasmo.add_edge(g, n2, n3)
		@linkconstraint(g, [i in Investments], n2[:AmountInvestedCurrentStage][i] == n3[:AmountInvestedPreviousStage][i])
	end
end

model = create_jump_graph_model(g)
model.solver = CplexSolver()
solve(model)
PlasmoAlgorithms.crosssolve(g; max_iterations=50, subgradientiterations=5, ϵ=0.005, timelimit=200, bdcuts=[:LP])
# bendersolve(g; max_iterations=50)






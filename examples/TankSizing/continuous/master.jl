function generate_master()
	m = Model(solver=BaronSolver(maxtime=5e4, epsr= 1e-3, prlevel=0, CplexLibName = "/opt/ibm/ILOG/CPLEX_Studio127/cplex/bin/x86-64_linux/libcplex1270.so"))
	# m = Model(solver=IpoptSolver())
	@variable(m, InventoryLowerBound[p]<=productTankSize[p in products]<=InventoryUpperBound[p])
	@variable(m, firstStageCost)

	@NLconstraint(m,  firstStageCost >= VariableInvestmentCostFactor * sum(prob[h] /TotalDemandPerDay[h] for h in scenarios) *sum(productTankSize[p]^0.5 for p in products))
	
	@objective(m, Min, firstStageCost )
	return m 

end























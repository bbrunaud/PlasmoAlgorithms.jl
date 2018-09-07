using Ipopt
include("input.jl")
function generate_model()
	m = Model(solver=BaronSolver(maxtime=5e4, epsr= 1e-3, CplexLibName = "/opt/ibm/ILOG/CPLEX_Studio127/cplex/bin/x86-64_linux/libcplex1270.so"))
	# m = Model(solver=IpoptSolver())
	@variable(m, 0<=campaignDuration[n in events, h in scenarios]<= maximum(ProductionLength_upper))
	@variable(m, 0<=amtProductInCampaign[p in products, n in events, h in scenarios]<=maximum(ProductionLength_upper.*MaxProductionRate))
	@variable(m, InventoryLowerBound[p]<=productInventory[p in products, n in events, h in scenarios]<=InventoryUpperBound[p])
	productInventory = getindex(m, :productInventory)
	setupperbound.(productInventory[1,1,:], InventoryLowerBound[1])
	setlowerbound.(productInventory[1,1,:], InventoryLowerBound[1])
	@variable(m, InventoryLowerBound[p]<=productTankSize[p in products]<=InventoryUpperBound[p])
	@variable(m, 0<=auxiliaryVariable[p in products, n in events, h in scenarios]<=InventoryUpperBound[p] - InventoryLowerBound[p])
	@variable(m, investmentCost[h in scenarios])
	@variable(m, setupCost[h in scenarios])
	@variable(m, variableCost[h in scenarios])
	@variable(m, 0<=cycleTime[h in scenarios]<=length(events)*maximum(ProductionLength_upper)+sum(CampaignSetupTime))
	@variable(m, costPerTon[h in scenarios])
	@variable(m, blin1[p in products, n in events, h in scenarios])

	@variable(m, assignProductToCampaign[p in products, n in events, h in scenarios], Bin)
	#the inital storage has some implications
	assignProductToCampaign = getindex(m, :assignProductToCampaign)
	setupperbound.(assignProductToCampaign[:, 1, :], 0)
	setlowerbound.(assignProductToCampaign[:, 1, :], 0)
	setupperbound.(assignProductToCampaign[1,1,:],  1)
	setlowerbound.(assignProductToCampaign[1,1,:], 1)
	setupperbound.(assignProductToCampaign[1,2,:], 0)
	setlowerbound.(assignProductToCampaign[1,2,:], 0)

	#constraints
	@constraint(m, TimeCapacity[h in scenarios], cycleTime[h] == sum( campaignDuration[n,h] + sum(CampaignSetupTime[p] * assignProductToCampaign[p,n,h] for p in products) for n in events))
	@constraint(m, UniqueProduct[n in events, h in scenarios], sum(assignProductToCampaign[p,n,h] for p in products) <= 1)

	@constraint(m, MaterialBalance[p in products, n in 1:(length(events)-1), h in scenarios], productInventory[p,n+1,h] == productInventory[p,n,h] + amtProductInCampaign[p,n,h] - DemandPerDay[p,h] * (campaignDuration[n,h] + sum(CampaignSetupTime[pp] * assignProductToCampaign[pp,n,h] for pp in products)))
	n = length(events)
	@constraint(m, MaterialBalance0[p in products, n, h in scenarios], productInventory[p,1,h] == productInventory[p,n,h] + amtProductInCampaign[p,n,h] - DemandPerDay[p,h] * (campaignDuration[n,h] + sum(CampaignSetupTime[pp] * assignProductToCampaign[pp,n,h] for pp in products)))


	@constraint(m, TankCapacity[p in products, n in events, h in scenarios], productInventory[p,n,h] <= productTankSize[p])
	@NLconstraint(m, blincon1[p in products, n in events, h in scenarios], blin1[p,n,h] == campaignDuration[n,h] * assignProductToCampaign[p,n,h] )
	@constraint(m, ProductionUpperBound[p in products, n in events, h in scenarios], amtProductInCampaign[p,n,h] <= MaxProductionRate[p] * blin1[p,n,h])
	@constraint(m, ProductionLowerBound[p in products, n in events, h in scenarios], amtProductInCampaign[p,n,h] >= MinProductionRate[p] * blin1[p,n,h])
	@constraint(m, CampaignUpperBound[n in events, h in scenarios], campaignDuration[n,h] <= sum(ProductionLength_upper[p] * assignProductToCampaign[p,n,h] for p in products))
	@constraint(m, CampaignLowerBound[n in events, h in scenarios], campaignDuration[n,h] >= sum(ProductionLength_lower[p] * assignProductToCampaign[p,n,h] for p in products))
	@constraint(m, CampaignSetupCostCon[h in scenarios], setupCost[h] == sum(CampaignSetupCost[p] * assignProductToCampaign[p,n,h] for p in products for n in events))
	@NLconstraint(m, CampaignInvestmentCost[h in scenarios], investmentCost[h] == VariableInvestmentCostFactor * sum(productTankSize[p]^0.5 for p in products))
	@NLconstraint(m, CampaignStorageCost[h in scenarios], variableCost[h] == sum(CampaignVariableCost[p] * auxiliaryVariable[p,n,h] * (campaignDuration[n,h] + sum(CampaignSetupTime[pp] * assignProductToCampaign[pp,n,h] for pp in products)) for p in products for n in events ))

	@constraint(m, AuxiliaryCon[p in products, n in 1:(length(events)-1), h in scenarios], auxiliaryVariable[p,n,h] == 0.5 * (productInventory[p,n+1,h] + productInventory[p,n,h]) - InventoryLowerBound[p])
	n = length(events)
	@constraint(m, AuxiliaryCon0[p in products, n, h in scenarios], auxiliaryVariable[p,n,h] == 0.5 * (productInventory[p,1,h] + productInventory[p,n,h]) - InventoryLowerBound[p])


	@constraint(m, CampaignCostPerTon[h in scenarios], costPerTon[h]*cycleTime[h]*TotalDemandPerDay[h] == investmentCost[h]*cycleTime[h] + setupCost[h] + variableCost[h])
	@constraint(m, Sequence[p in products, n in 1:(length(events)-1), h in scenarios], 1 - assignProductToCampaign[p,n,h] >= assignProductToCampaign[p,n+1,h])

	@constraint(m, BreakSymmetry[n in 1:(length(events)-1), h in scenarios], sum(assignProductToCampaign[p,n,h] for p in products) >= sum(assignProductToCampaign[p,n+1,h] for p in products))


	@objective(m, Min, sum(prob[h] * costPerTon[h] for h in scenarios))

	return m 

end

model = generate_model()
solve(model)






















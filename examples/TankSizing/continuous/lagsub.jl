function generate_lagsub(;prob = prob[1], DemandPerDay = DemandPerDay[:, 1], TotalDemandPerDay=TotalDemandPerDay[1], costPerTon_ub = 200, costPerTon_lb = 0, cycleTime_ub = 200, cycleTime_lb = 0, all_prob = prob, all_TotalDemandPerDay = TotalDemandPerDay )
	m = Model(solver=BaronSolver(maxtime=5e4, prlevel=0, epsr= 1e-4, CplexLibName = "/opt/ibm/ILOG/CPLEX_Studio127/cplex/bin/x86-64_linux/libcplex1270.so"))
	# m = Model(solver=IpoptSolver())
	@variable(m, InventoryLowerBound[p]<=productTankSize[p in products]<=InventoryUpperBound[p])
	@variable(m, firstStageCost)
	@variable(m, 0<=campaignDuration[n in events]<= maximum(ProductionLength_upper))
	@variable(m, 0<=amtProductInCampaign[p in products, n in events]<= ProductionLength_upper[p] * MaxProductionRate[p])
	@variable(m, InventoryLowerBound[p]<=productInventory[p in products, n in events]<=InventoryUpperBound[p])
	productInventory = getindex(m, :productInventory)
	setupperbound.(productInventory[1,1], InventoryLowerBound[1])
	setlowerbound.(productInventory[1,1], InventoryLowerBound[1])
	@variable(m, 0<=auxiliaryVariable[p in products, n in events]<=(InventoryUpperBound[p] - InventoryLowerBound[p])/30)
	@variable(m, investmentCost)
	@variable(m, setupCost)
	@variable(m, variableCost)
	@variable(m, cycleTime_lb<=cycleTime<=cycleTime_ub)
	@variable(m, costPerTon_lb<=costPerTon<=costPerTon_ub)
	#campaignLength is campaignDuration + setuptime 
	@variable(m, minimum(ProductionLength_lower + CampaignSetupTime)<=campaignLength[n in events]<= maximum(ProductionLength_upper + CampaignSetupTime)) 


	@variable(m, assignProductToCampaign[p in products, n in events], Bin)
	#the inital storage has some implications
	assignProductToCampaign = getindex(m, :assignProductToCampaign)
	setupperbound.(assignProductToCampaign[:, 1], 0)
	setlowerbound.(assignProductToCampaign[:, 1], 0)
	setupperbound.(assignProductToCampaign[1,1],  1)
	setlowerbound.(assignProductToCampaign[1,1], 1)
	setupperbound.(assignProductToCampaign[1,2], 0)
	setlowerbound.(assignProductToCampaign[1,2], 0)

	#add slack variable to make the problem have relatively complete recourse 
	@variable(m, slack[p in products]>=0)
	#constraints
	@constraint(m, TimeCapacity, cycleTime == sum( campaignDuration[n] + sum(CampaignSetupTime[p] * assignProductToCampaign[p,n] for p in products) for n in events))
	@constraint(m, UniqueProduct[n in events], sum(assignProductToCampaign[p,n] for p in products) <= 1)

	@constraint(m, MaterialBalance[p in products, n in 1:(length(events)-1)], productInventory[p,n+1] == productInventory[p,n] + amtProductInCampaign[p,n] - DemandPerDay[p] * (campaignDuration[n] + sum(CampaignSetupTime[pp] * assignProductToCampaign[pp,n] for pp in products)))
	n = length(events)
	@constraint(m, MaterialBalance0[p in products, n], productInventory[p,1] == productInventory[p,n] + amtProductInCampaign[p,n] - DemandPerDay[p] * (campaignDuration[n] + sum(CampaignSetupTime[pp] * assignProductToCampaign[pp,n] for pp in products)))


	@constraint(m, TankCapacity[p in products, n in events], productInventory[p,n] <= productTankSize[p] + slack[p])
	@constraint(m, ProductionUpperBound1[p in products, n in events], amtProductInCampaign[p,n] <= MaxProductionRate[p] * ProductionLength_upper[p] * assignProductToCampaign[p,n])
	@constraint(m, ProductionLowerBound1[p in products, n in events], amtProductInCampaign[p,n] >= MinProductionRate[p] * ProductionLength_lower[p] * assignProductToCampaign[p,n])
	@constraint(m, ProductionUpperBound2[p in products, n in events], amtProductInCampaign[p,n] <= MaxProductionRate[p] * campaignDuration[n])
	@constraint(m, ProductionLowerBound2[p in products, n in events], amtProductInCampaign[p,n] >= MinProductionRate[p] * (campaignDuration[n] - ProductionLength_upper[p] * (1 - assignProductToCampaign[p,n])))

	@constraint(m, CampaignUpperBound[n in events], campaignDuration[n] <= sum(ProductionLength_upper[p] * assignProductToCampaign[p,n] for p in products))
	@constraint(m, CampaignLowerBound[n in events], campaignDuration[n] >= sum(ProductionLength_lower[p] * assignProductToCampaign[p,n] for p in products))
	@constraint(m, CampanLengthDef[n in events], campaignLength[n] == campaignDuration[n] + sum(CampaignSetupTime[pp] * assignProductToCampaign[pp,n] for pp in products))
	@constraint(m, CampaignSetupCostCon, setupCost == sum(CampaignSetupCost[p] * assignProductToCampaign[p,n] for p in products for n in events))
	#@NLconstraint(m, CampaignInvestmentCost, investmentCost == VariableInvestmentCostFactor * sum(productTankSize[p]^0.5 for p in products))
	@NLconstraint(m, CampaignStorageCost, variableCost == sum(CampaignVariableCost[p] * auxiliaryVariable[p,n] * campaignLength[n] for p in products for n in events ))

	@constraint(m, AuxiliaryCon[p in products, n in 1:(length(events)-1)], auxiliaryVariable[p,n] == 0.5 * (productInventory[p,n+1] + productInventory[p,n]) - InventoryLowerBound[p])
	n = length(events)
	@constraint(m, AuxiliaryCon0[p in products, n], auxiliaryVariable[p,n] == 0.5 * (productInventory[p,1] + productInventory[p,n]) - InventoryLowerBound[p])

	#to do 
	@NLconstraint(m, CampaignCostPerTon, costPerTon*cycleTime*TotalDemandPerDay == setupCost + variableCost)
	@constraint(m, Sequence[p in products, n in 1:(length(events)-1)], 1 - assignProductToCampaign[p,n] >= assignProductToCampaign[p,n+1])

	@constraint(m, BreakSymmetry[n in 1:(length(events)-1)], sum(assignProductToCampaign[p,n] for p in products) >= sum(assignProductToCampaign[p,n+1] for p in products))


	@NLconstraint(m,  firstStageCost >= VariableInvestmentCostFactor * sum(all_prob[h] /all_TotalDemandPerDay[h] for h in scenarios) *sum(productTankSize[p]^0.5 for p in products))
	@objective(m, Min, prob * costPerTon + prob * sum(slack[p] for p in products) + prob * firstStageCost )

	return m 

end























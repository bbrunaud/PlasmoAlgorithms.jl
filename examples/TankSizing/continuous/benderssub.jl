function generate_benderssub(;prob = prob[1], DemandPerDay = DemandPerDay[:, 1], TotalDemandPerDay=TotalDemandPerDay[1], costPerTon_ub = 200, costPerTon_lb = 0, cycleTime_ub = 200, cycleTime_lb = 0 )
	m = Model(solver=CplexSolver(CPX_PARAM_SCRIND=0))
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

	#bilinear terms that represent auxiliaryVariable[p,n] * campaignLength[n]
	@variable(m, bilinear[p in products, n in events]>=0)
	#bilinear term that represent costPerTon * cycleTime
	@variable(m, costMultCycleTime)

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
	@constraint(m, CampaignStorageCost, variableCost == sum(CampaignVariableCost[p] * bilinear[p,n] for p in products for n in events ))
	# ========piecewise Mccormick envelopes
	#auxiliaryVariable[p,n] * campaignLength[n]
	points = 1:2
	intervals = 1:(length(points)-1)
	auxiliaryVariable_lb = zeros(length(products), length(events))
	auxiliaryVariable_ub = zeros(length(products), length(events))
	for n in events
		auxiliaryVariable_ub[:, n] = (InventoryUpperBound - InventoryLowerBound)/30
	end
	campaignLength_lb =  minimum(ProductionLength_lower + CampaignSetupTime)
	campaignLength_ub =  maximum(ProductionLength_upper + CampaignSetupTime)
	campaignLength_points = ones(length(points)) * campaignLength_lb
	for t in intervals
		campaignLength_points[t+1] = campaignLength_lb + t * (campaignLength_ub - campaignLength_lb) / length(intervals)
	end
	@variable(m, dot_auxiliaryVariable[p in products, n in events, t in intervals]>=0)
	@variable(m, dot_bilinear[p in products, n in events, t in intervals]>=0)
	@variable(m, dot_campaignLength[n in events, t in intervals]>=0)
	# @variable(m, beta[n in events, t in intervals], Bin)
	@variable(m, 0<=beta[n in events, t in intervals]<=1)

	@constraint(m, n1[p in products, n in events], sum(dot_auxiliaryVariable[p,n,t] for t in intervals) == auxiliaryVariable[p,n])
	@constraint(m, n2[p in products, n in events], sum(dot_bilinear[p,n,t] for t in intervals) == bilinear[p,n])
	@constraint(m, n3[n in events], sum(dot_campaignLength[n,t] for t in intervals) == campaignLength[n])
	@constraint(m, n4[n in events], sum(beta[n,t] for t in intervals) == 1)
	
	@constraint(m, Mc1[p in products, n in events, t in intervals], dot_bilinear[p,n,t] >= auxiliaryVariable_lb[p,n] * dot_campaignLength[n,t] + dot_auxiliaryVariable[p,n,t] * campaignLength_points[t] - auxiliaryVariable_lb[p,n] * campaignLength_points[t] * beta[n,t])
	@constraint(m, Mc2[p in products, n in events, t in intervals], dot_bilinear[p,n,t] >= auxiliaryVariable_ub[p,n] * dot_campaignLength[n,t] + dot_auxiliaryVariable[p,n,t] * campaignLength_points[t+1] - auxiliaryVariable_ub[p,n] * campaignLength_points[t+1] * beta[n,t])
	@constraint(m, Mc3[p in products, n in events, t in intervals], dot_bilinear[p,n,t] <= auxiliaryVariable_ub[p,n] * dot_campaignLength[n,t] + dot_auxiliaryVariable[p,n,t] * campaignLength_points[t] - auxiliaryVariable_ub[p,n] * campaignLength_points[t] * beta[n,t])
	@constraint(m, Mc4[p in products, n in events, t in intervals], dot_bilinear[p,n,t] <= auxiliaryVariable_lb[p,n] * dot_campaignLength[n,t] + dot_auxiliaryVariable[p,n,t] * campaignLength_points[t+1] - auxiliaryVariable_lb[p,n] * campaignLength_points[t+1] * beta[n,t])
	@constraint(m, b4[n in events, t in intervals], dot_campaignLength[n,t] <= campaignLength_points[t+1] * beta[n,t])
	@constraint(m, b5[n in events, t in intervals], dot_campaignLength[n,t] >= campaignLength_points[t] * beta[n,t])
	@constraint(m, b6[p in products, n in events, t in intervals], dot_auxiliaryVariable[p,n,t] <= beta[n,t] * auxiliaryVariable_ub[p,n])

	# =======end 

	@constraint(m, AuxiliaryCon[p in products, n in 1:(length(events)-1)], auxiliaryVariable[p,n] == 0.5 * (productInventory[p,n+1] + productInventory[p,n]) - InventoryLowerBound[p])
	n = length(events)
	@constraint(m, AuxiliaryCon0[p in products, n], auxiliaryVariable[p,n] == 0.5 * (productInventory[p,1] + productInventory[p,n]) - InventoryLowerBound[p])

	#to do 
	@constraint(m, CampaignCostPerTon, costMultCycleTime*TotalDemandPerDay == setupCost + variableCost)
	#piecewise Mccormick Envelope for costMultCycleTime=======
	#discretize cycleTime
	cycleTime_points = ones(length(points)) * cycleTime_lb
	for t in 1:(length(points)-1)
		cycleTime_points[t+1] = cycleTime_lb + t * (cycleTime_ub - cycleTime_lb) / (length(points) - 1)
	end

	# @variable(m, delta[t in intervals], Bin)
	@variable(m, 0<=delta[t in intervals]<=1)
	@variable(m, dot_cycleTime[t in intervals]>=0)
	@variable(m, dot_costPerTon[t in intervals]>=0)
	@variable(m, dot_costMultCycleTime[t in intervals]>=0)
	@constraint(m, sum(delta[t] for t in intervals) == 1)
	@constraint(m, sum(dot_cycleTime[t] for t in intervals) == cycleTime)
	@constraint(m, sum(dot_costPerTon[t] for t in intervals) == costPerTon)
	@constraint(m, sum(dot_costMultCycleTime[t] for t in intervals) == costMultCycleTime)

	@constraint(m, mc1[t in intervals], dot_costMultCycleTime[t] >= costPerTon_lb * dot_cycleTime[t] + dot_costPerTon[t] * cycleTime_points[t] - costPerTon_lb * cycleTime_points[t] * delta[t])
	@constraint(m, mc2[t in intervals], dot_costMultCycleTime[t] >= costPerTon_ub * dot_cycleTime[t] + dot_costPerTon[t] * cycleTime_points[t+1] - costPerTon_ub * cycleTime_points[t+1] * delta[t])
	@constraint(m, mc3[t in intervals], dot_costMultCycleTime[t] <= costPerTon_ub * dot_cycleTime[t] + dot_costPerTon[t] * cycleTime_points[t] - costPerTon_ub * cycleTime_points[t] * delta[t])
	@constraint(m, mc4[t in intervals], dot_costMultCycleTime[t] <= costPerTon_lb * dot_cycleTime[t] + dot_costPerTon[t] * cycleTime_points[t+1] - costPerTon_lb * cycleTime_points[t+1] * delta[t])
	@constraint(m, b1[t in intervals], dot_cycleTime[t] <= cycleTime_points[t+1]*delta[t])
	@constraint(m, b2[t in intervals], dot_cycleTime[t] >= cycleTime_points[t]*delta[t])
	@constraint(m, b3[t in intervals], dot_costPerTon[t] <= costPerTon_ub * delta[t])
	# ===========end

	@constraint(m, Sequence[p in products, n in 1:(length(events)-1)], 1 - assignProductToCampaign[p,n] >= assignProductToCampaign[p,n+1])

	@constraint(m, BreakSymmetry[n in 1:(length(events)-1)], sum(assignProductToCampaign[p,n] for p in products) >= sum(assignProductToCampaign[p,n+1] for p in products))


	@objective(m, Min, prob * costPerTon + prob * sum(slack[p] for p in products)  )

	return m 

end























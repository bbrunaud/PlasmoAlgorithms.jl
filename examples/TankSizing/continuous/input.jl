using SpecialFunctions
using JuMP
using BARON

#sets
products = 1:3
events = 1:3
scenarios = 1:3
realizations = 1:3
# NOTE: Use |scenarios| = |realizations| or |realizations|^2 or |realizations|^3


# *-------------------------------------------------
# *		DEFINE MODEL-SPECIFIC PARAMETERS
# *-------------------------------------------------
VariableInvestmentCostFactor = 0.3271
NumDaysInYear = 365

MinProductionRate = 
	[	15.0
			15.0
			7.0]

MaxProductionRate =
	[50.0
	50.0
	50.0]

InventoryLowerBound =
	[643.0
	536.0
	214.0]

InventoryUpperBound =
	[4018.36
	3348.63
	1339.45]

InitialInventory = zeros(length(products))

ProductionLength_lower = 
	[1
	1
	1]

ProductionLength_upper =
	[40
	40
	40]

ProductDemand_nominal = 
	[4190
	3492
	1397]

ProductDemand_stdev = zeros(length(products))

ProductDemand = zeros(length(products),length(scenarios)) 

CampaignSetupTime = 
	[0.4
	0.2
	0.1]

CampaignVariableCost = 
	[18.8304
	19.2934
	19.7563 ]

CampaignSetupCost =
	[10
	20
	30]

prob = zeros(length(scenarios))
subprob = zeros(length(realizations)) 

# *---------------------------------------------
# *	 Compute derived parameter values
# *---------------------------------------------

InitialInventory = 1.1*InventoryLowerBound
ProductDemand_stdev = 0.1 * ProductDemand_nominal


# =========== Generate scenarios for the uncertain parameters =============
#initialize product demand 
for index in scenarios 
	ProductDemand[:, index] = ProductDemand_nominal[:]
end

function errorf(x)
	return (1 + erf(x/sqrt(2)))/2
end

for index in realizations
	if length(realizations) > 1 && index == 1
		subprob[index] = errorf(-3+6/length(realizations))
	end
	if length(realizations) > 1 && index>1 && index < length(realizations)
		subprob[index] = errorf(-3+index*6/length(realizations))-errorf(-3+(index-1)*6/length(realizations))
	end 
	if length(realizations) > 1 && index == length(realizations)
		subprob[index] = 1 - errorf(-3+(length(realizations)-1)*6/length(realizations))
	end
end

if length(scenarios) == 1
	prob[1] = 1
elseif length(scenarios) == length(realizations)
	for index in realizations
		prob[index] = subprob[index]
		ProductDemand[1, index] = ProductDemand_nominal[1] - 3 * ProductDemand_stdev[1] * (1 - 1/length(realizations))  + (index - 1) * 6 * ProductDemand_stdev[1] / length(realizations)
	end
elseif length(scenarios) == length(realizations) * length(realizations)
	for index1 in realizations
		for index2 in realizations
			index = index1 + (index2 - 1) * length(realizations)
			prob[index] = subprob[index1] * subprob[index2]
			ProductDemand[1, index] = ProductDemand_nominal[1] - 3 * ProductDemand_stdev[1] * (1 - 1/length(realizations))  + (index1 - 1) * 6 * ProductDemand_stdev[1] / length(realizations)
			ProductDemand[2, index] = ProductDemand_nominal[2] - 3 * ProductDemand_stdev[2] * (1 - 1/length(realizations))  + (index2 - 1) * 6 * ProductDemand_stdev[2] / length(realizations)
		end
	end
elseif length(scenarios) == length(realizations) * length(realizations) * length(realizations)
	for index1 in realizations
		for index2 in realizations
			for index3 in realizations
				index = index1 + (index2 - 1) * length(realizations) + (index3 -1 ) * length(realizations) * length(realizations)
				prob[index] = subprob[index1] * subprob[index2] * subprob[index3]
				ProductDemand[1, index] = ProductDemand_nominal[1] - 3 * ProductDemand_stdev[1] * (1 - 1/length(realizations))  + (index1 - 1) * 6 * ProductDemand_stdev[1] / length(realizations)
				ProductDemand[2, index] = ProductDemand_nominal[2] - 3 * ProductDemand_stdev[2] * (1 - 1/length(realizations))  + (index2 - 1) * 6 * ProductDemand_stdev[2] / length(realizations)
				ProductDemand[3, index] = ProductDemand_nominal[3] - 3 * ProductDemand_stdev[3] * (1 - 1/length(realizations))  + (index3 - 1) * 6 * ProductDemand_stdev[3] / length(realizations)
			end
		end
	end	
else 
	error("ERROR in setting the number of scenarios! Try again")
end

# *---------------------------------------------
# *	 Compute derived parameter values
# *---------------------------------------------
DemandPerDay = zeros(length(products), length(scenarios))
TotalDemandPerDay = zeros(length(scenarios))
DemandPerDay = ProductDemand/NumDaysInYear
TotalDemandPerDay[:] = sum(DemandPerDay[i, :] for i in 1:length(products));
CampaignVariableCost = CampaignVariableCost/NumDaysInYear;




























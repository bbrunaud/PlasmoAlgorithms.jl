$ONTEXT

Tank sizing problem based on the three product example in 
Rebennack et al., Computers and Chemical Engineering, 2011
This model contains three uncertain product demands.
Number of binary complicating variables: 
Number of continuous complicating variables: 
Number of binary recourse variables: 0
Number of continuous recourse variables: *s
Number of bilinear terms: *s
Number of univariate signomial terms: 
Number of complicating constraints: 
Number of recourse constraints: 

This file is based on tanksize.350 in the gamslib_ml folder of the GAMS
distribution

$OFFTEXT

*-------------------------------------------
*			SET SOLUTION OPTIONS
*-------------------------------------------

*OPTION LIMROW = 0;
*OPTION LIMCOL = 0;
OPTION OPTCA  = 1E-09;
OPTION OPTCR  = 1E-03;
OPTION RESLIM = 3E+02;
OPTION ITERLIM = 1E+09;

OPTION LP=CPLEX;
OPTION NLP=SNOPT;
OPTION MIP=CPLEX;
*OPTION MINLP=bonmin;
*OPTION MINLP=ANTIGONE;
OPTION MINLP=BARON;
*OPTION MINLP=COUENNE;
*OPTION MINLP=SCIP;

*--------------------------------
*		SET DEFINITIONS
*--------------------------------

SETS
	p		"products"									/ 1*3 /
	n		"event points"								/ 1*3 /
	h		"number of scenarios"						/ 1*1 /	
	subh	"num. realizations per uncertain param"		/ 1*1 /	
;
*** NOTE: Use h = subh or subh^2 or subh^3

alias(subh,subh2,subh3);
alias(p,pp);

*-------------------------------------------------
*		DEFINE MODEL-SPECIFIC PARAMETERS
*-------------------------------------------------

SCALARS
	VariableInvestmentCostFactor	"variable part of the tank investment cost"	/ 0.3271 /
	NumDaysInYear "number of days in a year" / 365 /
;

PARAMETERS
	MinProductionRate(p) "lower bound on the production rate in m^3/day"
	/	1	15.0
		2	15.0
		3	7.0/,

	MaxProductionRate(p) "upper bound on the production rate in m^3/day"
	/	1	50.0
		2	50.0
		3	50.0/,

	InventoryLowerBound(p) "lower bound on inventory in m^3"
	/	1	643.0
		2	536.0
		3	214.0/,

	InventoryUpperBound(p) "upper bound on inventory in m^3"
	/	1	4018.36
		2	3348.63
		3	1339.45/,

	InitialInventory(p) "initial inventory in m^3",

	ProductionLength_lower(p) "lower bound on production length"
	/	1	1
		2	1
		3	1/,

	ProductionLength_upper(p) "upper bound on production length"
	/	1	40
		2	40
		3	40/,

	ProductDemand_nominal(p) "nominal demand of product in m^3/year"
	/	1	4190
		2	3492
		3	1397/,

	ProductDemand_stdev(p) "standard deviation of demand of product in m^3/year",

	ProductDemand(p,h) "demand of product in m^3/year",

	CampaignSetupTime(p) "campaign setup time in days"
	/	1	0.4
		2	0.2
		3	0.1/,

	CampaignVariableCost(p) "tank variable cost per ton"
	/	1	18.8304
		2	19.2934
		3	19.7563/,

	CampaignSetupCost(p) "campaign setup cost"
	/	1	10
		2	20
		3	30/,

	prob(h) "probability of each scenario",

	subprob(subh) "probability for each individual uncertain realization"
;

*---------------------------------------------
*	 Compute derived parameter values
*---------------------------------------------

InitialInventory(p) = 1.1*InventoryLowerBound(p);

ProductDemand_stdev('1') = ProductDemand_nominal('1')*0.1;
ProductDemand_stdev('2') = ProductDemand_nominal('2')*0.1;
ProductDemand_stdev('3') = ProductDemand_nominal('3')*0.1;

*---------------------------------------------
*	 Initialize the uncertain parameters
*---------------------------------------------

ProductDemand(p,h) = ProductDemand_nominal(p);



*=========== Generate scenarios for the uncertain parameters =============


if (card(h)=1,
	prob(h)=1;

elseif (card(h)=card(subh)),
	subprob(subh)$(card(subh)>1 and ord(subh)=1) = errorf(-3+6/card(subh));
	subprob(subh)$(card(subh)>1 and ord(subh)>1 and ord(subh)<card(subh))
		= errorf(-3+ord(subh)*6/card(subh))-errorf(-3+(ord(subh)-1)*6/card(subh));
	subprob(subh)$(card(subh)>1 and ord(subh) = card(subh))
		= 1 - errorf(-3+(card(subh)-1)*6/card(subh));

	loop(subh,
		prob(h)$(ord(h)=ord(subh)) = subprob(subh);

		ProductDemand('1',h)$(ord(h)=ord(subh)) = 	ProductDemand_nominal('1') -
													3*ProductDemand_stdev('1')*(1 - 1/card(subh)) +
													(ord(subh)-1)*6*ProductDemand_stdev('1')/card(subh);
	);

elseif (card(h)=card(subh)*card(subh)),
	subprob(subh)$(card(subh)>1 and ord(subh)=1) = errorf(-3+6/card(subh));
	subprob(subh)$(card(subh)>1 and ord(subh)>1 and ord(subh)<card(subh))
		= errorf(-3+ord(subh)*6/card(subh))-errorf(-3+(ord(subh)-1)*6/card(subh));
	subprob(subh)$(card(subh)>1 and ord(subh) = card(subh))
		= 1 - errorf(-3+(card(subh)-1)*6/card(subh));

	loop(subh,
	    loop(subh2,
			prob(h)$(ord(h)=ord(subh)+(ord(subh2)-1)*card(subh))
			= subprob(subh)*subprob(subh2);

			ProductDemand('1',h)$(ord(h)=ord(subh)+(ord(subh2)-1)*card(subh))
			= ProductDemand_nominal('1') - 3*ProductDemand_stdev('1')*(1 - 1/card(subh)) +
				(ord(subh)-1)*6*ProductDemand_stdev('1')/card(subh);

			ProductDemand('2',h)$(ord(h)=ord(subh)+(ord(subh2)-1)*card(subh))
			= ProductDemand_nominal('2') - 3*ProductDemand_stdev('2')*(1 - 1/card(subh2)) +
				(ord(subh2)-1)*6*ProductDemand_stdev('2')/card(subh2);
	    );
	);

elseif (card(h)=card(subh)*card(subh)*card(subh)),
	subprob(subh)$(card(subh)>1 and ord(subh)=1) = errorf(-3+6/card(subh));
	subprob(subh)$(card(subh)>1 and ord(subh)>1 and ord(subh)<card(subh))
		= errorf(-3+ord(subh)*6/card(subh))-errorf(-3+(ord(subh)-1)*6/card(subh));
	subprob(subh)$(card(subh)>1 and ord(subh) = card(subh))
		= 1 - errorf(-3+(card(subh)-1)*6/card(subh));

	loop(subh,
	    loop(subh2,
			loop(subh3,
				prob(h)$(ord(h)=ord(subh)+(ord(subh2)-1)*card(subh)+(ord(subh3)-1)*card(subh)*card(subh2))
				= subprob(subh)*subprob(subh2)*subprob(subh3);

				ProductDemand('1',h)$(ord(h)=ord(subh)+(ord(subh2)-1)*card(subh)+(ord(subh3)-1)*card(subh)*card(subh2))
				= ProductDemand_nominal('1') - 3*ProductDemand_stdev('1')*(1 - 1/card(subh)) +
					(ord(subh)-1)*6*ProductDemand_stdev('1')/card(subh);

				ProductDemand('2',h)$(ord(h)=ord(subh)+(ord(subh2)-1)*card(subh)+(ord(subh3)-1)*card(subh)*card(subh2))
				= ProductDemand_nominal('2') - 3*ProductDemand_stdev('2')*(1 - 1/card(subh2)) +
					(ord(subh2)-1)*6*ProductDemand_stdev('2')/card(subh2);

				ProductDemand('3',h)$(ord(h)=ord(subh)+(ord(subh2)-1)*card(subh)+(ord(subh3)-1)*card(subh)*card(subh2))
				= ProductDemand_nominal('3') - 3*ProductDemand_stdev('3')*(1 - 1/card(subh3)) +
					(ord(subh3)-1)*6*ProductDemand_stdev('3')/card(subh3);
			);
	    );
	);

else
	abort "ERROR in setting the number of scenarios! Try again"
);

*======================================================================

*---------------------------------------------
*	 Compute derived parameter values
*---------------------------------------------

PARAMETERS
	DemandPerDay(p,h)		"demand/day/product [tons/day]"
	TotalDemandPerDay(h)	"total demand/day [tons/day]"
;

DemandPerDay(p,h) = ProductDemand(p,h)/NumDaysInYear;
TotalDemandPerDay(h) = sum(p, DemandPerDay(p,h));
CampaignVariableCost(p) = CampaignVariableCost(p)/NumDaysInYear;


*---------------------------------------------------------------
*-------------    FULL SPACE PROBLEM DEFINITION   --------------
*---------------------------------------------------------------


POSITIVE VARIABLES
	campaignDuration(n,h)		"duration of the campaigns"
	amtProductInCampaign(p,n,h)	"amount of product p produced in campaign n"
	productInventory(p,n,h)		"amount of product p stored at the beginning of campaign n"
	productTankSize(p)			"size of the product tanks in tons"
	auxiliaryVariable(p,n,h)	"auxiliary variables"
	investmentCost(h)			"investment costs"
	setupCost(h)				"campaign setup costs"
	variableCost(h)				"variable storage costs"
	cycleTime(h)				"cycle time"
	costPerTon(h)				"cost per ton"
	blin1(p,n,h)				"bilinear terms"
;

BINARY VARIABLES
	assignProductToCampaign(p,n,h)	"binary variable mapping product to campaign"
;

VARIABLE
	objvar		"objective function"
;


EQUATIONS
	TimeCapacity(h)					"time capacity"
	UniqueProduct(n,h)				"at most one product per campaign"
	MaterialBalance(p,n,h)			"material balance constraint (steady state)"
	TankCapacity(p,n,h)				"tank capacity constraint"
	blincon1						"constraint defining bilinear terms"
	ProductionUpperBound(p,n,h)		"upper bound for product"
	ProductionLowerBound(p,n,h)		"lower bound for product"
	CampaignUpperBound(n,h)			"upper bound on duration"
	CampaignLowerBound(n,h)			"lower bound on duration"
	CampaignSetupCostCon(h)			"campaign setup cost"
	CampaignInvestmentCost(h)		"campaign investment cost"
	CampaignStorageCost(h)			"campaign variable storage cost"
	AuxiliaryCon(p,n,h)				"define the auxiliary variables"
	CampaignCostPerTon(h)			"cost per ton"
	Sequence(p,n,h)					"redundant constraint on the omega"
	BreakSymmetry(n,h)				"break the symmetry of active campaigns"
	objective						"objective function"
;

*-------------------------------------
*		EQUATION DEFINITIONS
*-------------------------------------

objective .. objvar =e= sum(h,prob(h)*costPerTon(h))
;

*** time balance constraint with unknown cycle time T
TimeCapacity(h) .. cycleTime(h) =e= sum(n, campaignDuration(n,h) + sum (p, CampaignSetupTime(p)*assignProductToCampaign(p,n,h))) 
;

UniqueProduct(n,h) .. sum(p, assignProductToCampaign(p,n,h)) =l= 1
;

MaterialBalance(p,n,h) .. productInventory(p,n++1,h) =e= productInventory(p,n,h) + amtProductInCampaign(p,n,h) - 
															DemandPerDay(p,h)*(campaignDuration(n,h) + 
																				sum (pp, CampaignSetupTime(pp)*assignProductToCampaign(pp,n,h)))
;

TankCapacity(p,n,h) .. productInventory(p,n,h) =l= productTankSize(p)
;

blincon1(p,n,h) .. blin1(p,n,h) =e= campaignDuration(n,h)*assignProductToCampaign(p,n,h)
;

ProductionUpperBound(p,n,h) .. amtProductInCampaign(p,n,h) =l= MaxProductionRate(p)*blin1(p,n,h)
;

ProductionLowerBound(p,n,h) .. amtProductInCampaign(p,n,h) =g= MinProductionRate(p)*blin1(p,n,h)
;

*** semi-continuous upper and lower bound on campaigns
CampaignUpperBound(n,h) .. campaignDuration(n,h) =l= sum(p, ProductionLength_upper(p)*assignProductToCampaign(p,n,h))
;

CampaignLowerBound(n,h) .. campaignDuration(n,h) =g= sum(p, ProductionLength_lower(p)*assignProductToCampaign(p,n,h))
;

CampaignSetupCostCon(h) .. setupCost(h) =e= sum((p,n), CampaignSetupCost(p)*assignProductToCampaign(p,n,h))
;

CampaignInvestmentCost(h) .. investmentCost(h) =e= VariableInvestmentCostFactor*sum(p, sqrt(productTankSize(p)))
;

CampaignStorageCost(h) .. variableCost(h) =e= sum((p,n), CampaignVariableCost(p)*auxiliaryVariable(p,n,h)*
													(campaignDuration(n,h) + sum(pp, CampaignSetupTime(pp)*assignProductToCampaign(pp,n,h))))
;

AuxiliaryCon(p,n,h) .. auxiliaryVariable(p,n,h) =e= 0.5*(productInventory(p,n++1,h) + productInventory(p,n,h)) - InventoryLowerBound(p)
;

CampaignCostPerTon(h) .. costPerTon(h)*cycleTime(h)*TotalDemandPerDay(h) =e= investmentCost(h)*cycleTime(h) + setupCost(h) + variableCost(h)
;

*** if a product is produced during period n, then it cannot be produced during period n+1
Sequence(p,n,h) .. 1 - assignProductToCampaign(p,n,h) =g= assignProductToCampaign(p,n+1,h)
;

*** break symmetry by shifting empty periods to the end
BreakSymmetry(n,h) .. sum(p, assignProductToCampaign(p,n,h)) =g= sum(p, assignProductToCampaign(p,n+1,h))
;


MODEL tanksize1 / all /;


*---------------------------------------------
*			SET VARIABLE BOUNDS
*---------------------------------------------

campaignDuration.up(n,h) = sum(p,ProductionLength_upper(p));

amtProductInCampaign.up(p,n,h) = MaxProductionRate(p)*campaignDuration.up(n,h);

productInventory.lo(p,n,h) = InventoryLowerBound(p);
productInventory.up(p,n,h) = InventoryUpperBound(p);
*** fix initial storage
productInventory.fx('1','1',h) = InventoryLowerBound('1');

auxiliaryVariable.up(p,n,h) = InventoryUpperBound(p) - InventoryLowerBound(p);

cycleTime.up(h) = card(n)*sum(p,ProductionLength_upper(p)) + sum(p,CampaignSetupTime(p));

*** the inital storage has some implications
assignProductToCampaign.fx(p,'1',h) = 0;
assignProductToCampaign.fx('1','1',h) = 1;
assignProductToCampaign.fx('1','2',h) = 0;

productTankSize.lo(p) = InventoryLowerBound(p);
productTankSize.up(p) = InventoryUpperBound(p);


*---------------------------------------------
*			SET VARIABLE BOUNDS
*---------------------------------------------

SOLVE tanksize1 minimizing objvar using MINLP;


*---------------------------------------------------
*			PRINT FINAL SOLUTION STATISTICS
*---------------------------------------------------

SCALARS SOLVER_TIME, WALL_TIME;

SOLVER_TIME = tanksize1.resusd;
WALL_TIME = tanksize1.etsolve;

display SOLVER_TIME, WALL_TIME;


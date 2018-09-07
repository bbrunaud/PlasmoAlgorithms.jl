*-------------------------------------------
*			SET SOLUTION OPTIONS
*-------------------------------------------

OPTION LIMROW = 0;
OPTION LIMCOL = 0;
OPTION OPTCA  = 1E-09;
OPTION OPTCR  = 1E-03;
OPTION RESLIM = 3E+02;
OPTION ITERLIM = 1E+09;

OPTION LP=CPLEX;
OPTION NLP=SNOPT;
OPTION MIP=CPLEX;
*OPTION MINLP=bonmin;
OPTION MINLP=ANTIGONE;
*OPTION MINLP=BARON;
*OPTION MINLP=COUENNE;
*OPTION MINLP=SCIP;

*--------------------------------
*		SET DEFINITIONS
*--------------------------------

SETS
	p		"products"									/ 1*3 /
	n		"event points"								/ 1*3 /
	h		"number of scenarios"						/ 1*2 /	
	subh	"num. realizations per uncertain param"		/ 1*2 /	
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
execute_unload "results_2.gdx" DemandPerDay TotalDemandPerDay;
*execute 'gdxxrw.exe results_1.gdx o=results_1.xls par=DemandPerDay TotalDemandPerDay'
execute 'gdx2xls results_2.gdx'

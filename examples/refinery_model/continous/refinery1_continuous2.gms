$ONTEXT

Bilinear refinery model based on the example in
Yang et al., AIChE J., 2016
This model contains uncertain Sulphur and VR yields.
Number of binary complicating variables: 
Number of continuous complicating variables: 0
Number of binary recourse variables: 0
Number of continuous recourse variables: *s
Number of bilinear terms: *s
Number of complicating constraints: 
Number of recourse constraints: 

$OFFTEXT

*--------------------------------
*		SET DEFINITIONS
*--------------------------------

SETS
	c		"crudes"								/ 1*10 /
	w		"components"							/ 1*8 /
	Re_in	"reformer input"						/ 1*2 /
	Re_out	"reformer output"						/ 1*6 /
	Cr_in	"cracker input"							/ 1*2 /
	Cr_out	"cracker output"						/ 1*6 /
	Cr_CGO	"cracker-CGO output"					/ 1*3 /
	Cr_mode	"cracker modes"							/ 1*2 /
	Iso_out	"isomerisation output"					/ 1*4 /
	Des_out	"desulphurisation output"				/ 1*4 /
	PG98_in	"PG98 input"							/ 1*6 /
	Burn	"burn streams"							/ 1*3 /
	JPF_in	"JPF input"								/ 1*3 /
	JPF_out	"JPF output"							/ 1*2 /
	AGO_in	"AGO input"								/ 1*3 /
	p		"products"								/ 1*7 /
	LG_in	"LG input"								/ 1*5 /
	LG_out	"LG output"								/ 1*4 /
	LG_prop	"LG properties"							/ 1*3 /
	h		"scenarios"								/ 1*5 /
;

alias(w,w2);

*-------------------------------------------
*			SET SOLUTION OPTIONS
*-------------------------------------------

OPTION LIMROW = 0;
OPTION LIMCOL = 0;
OPTION OPTCA  = 1E-09;
OPTION OPTCR  = 1E-05;
OPTION RESLIM = 1E+04;
OPTION ITERLIM = 1E+09;

OPTION LP=CPLEX;
*OPTION NLP=SNOPT;
OPTION MIP=CPLEX;
*OPTION MINLP=antigone;
OPTION MINLP=BARON;

*---------------------------------------------
*		DEFINE PROBLEM-SPECIFIC PARAMETERS
*---------------------------------------------

SCALARS
	Desulphurisation_capacity  / 125 /,
	CDU_capacity  / 700 /,
	Reformer95_lower  / 5 /,
	Reformer_capacity  / 65 /,
	Cracker_capacity  / 175 /,
	GranularityOfBarrels  / 5000 /,
	LG_sale  / 561.6 /,
	LN_sale  / 1003 /,
	HF_sale  / 637 /,
	ES95_sale  / 1194 /,
	PG98_sale  / 1231 /,
	JET_sale  / 923 /,
	AGO_sale  / 907 /,
	CGO_density  / 0.95 /,
	Mogas_viscosity  / 12.2 /,
	AGO_viscosity  / 11.65 /,
	Mogas_Sulphur  / 2.1 /,
	AGO_Sulphur  / 1.68 /,
	Desulphurisation_CGO_cost,
	Isomerisation_cost  / 6 /,
	Reformer95_cost  / 2.7 /,
	Reformer100_cost  / 3.2 /,
	Cracker_Mogas_cost  / 3.2 /,
	Cracker_AGO_cost  / 3 /,
	Barrel_lower_bound  / 100000 /,
	Barrel_upper_bound  / 1500000 /,
	Sulphur_spec  / 0.0015 /
;

Desulphurisation_CGO_cost = ((Mogas_Sulphur*109.0909 + 365.4546)/1000)*(0.85/0.159)/CGO_density;

TABLE Reformer_fraction(Re_in,Re_out)
		1     2     3     4     5     6
	1	0.08  0.09  0.83  0     0.019 2.7
	2	0.09  0.12  0     0.79  0.026 3.2
;

TABLE Cracker_fraction(Cr_in,Cr_out)
		1     2     3     4     5     6
	1	0.015 0.053 0.436 0.446 0.007 3.2
	2	0.012 0.046 0.381 0.511 0.007 3.0
;

TABLE Desulphurisation_fraction(c,Des_out)
		1     2     3     4
	1	0.98  0.02  0.02  0
	2	0.98  0.02  0.02  0
	3	0.98  0.02  0.02  0
	4	0.98  0.02  0.02  0
	5	0.98  0.02  0.02  0
	6	0.97  0.03  0.02  0
	7	0.97  0.03  0.02  0
	8	0.96  0.04  0.02  0
	9	0.98  0.02  0.02  0
	10	0.96  0.04  0.02  0
;

TABLE JPF_fraction(JPF_in,JPF_out)
		1     2
	1	0.05  0.035
	2	0.10  0.065
	3	0.85  0.900
;

TABLE Crude_yield(c,w)
		1      2      3      4      5      6      7      8
	1	0.0020 0.0091 0.0698 0.1598 0.1003 0.2876 0.2682 0.1032
	2	0.0020 0.0089 0.0480 0.0959 0.0796 0.2249 0.2735 0.2672
	3	0.0020 0.0080 0.0610 0.1206 0.0861 0.2414 0.2646 0.2163
	4	0.0040 0.0200 0.0851 0.1532 0.0947 0.2539 0.2535 0.1356
	5	0.0020 0.0115 0.0543 0.1026 0.0765 0.2286 0.2695 0.2550
	6	0.0010 0.0064 0.0246 0.0607 0.0518 0.1900 0.2932 0.3723
	7	0.0020 0.0155 0.0945 0.1661 0.1160 0.2656 0.2317 0.1086
	8	0.0029 0.0130 0.0652 0.1196 0.0838 0.2127 0.2408 0.2620
	9	0.0040 0.0157 0.0749 0.1267 0.0915 0.2353 0.2510 0.2009
	10	0.0040 0.0107 0.0604 0.1123 0.0784 0.2092 0.2491 0.2759
;

TABLE LG_parameters(LG_prop,LG_in)
		1     2     3     4     5
	1	4.30  4.28  4.36  4.21  4.22
	2	93.0  92.5  93.6  92.7  93.9
	3	90.0  89.6  90.9  89.0  90.2
;

*-------------------------------------------------
*		DEFINE PARAMETERS
*-------------------------------------------------

PARAMETERS

	Isomerisation_fraction(Iso_out)
	/	1	0.03
		2	0.97
		3	0.04
		4	6.0/,

	Desulphurisation_fraction2(Des_out)
	/	1	0.96
		2	0.04
		3	0.02
		4	20.0/,

	Crude_density(c)
	/	1	0.8441
		2	0.8910
		3	0.8441
		4	0.8369
		5	0.8829
		6	0.9315
		7	0.8252
		8	0.8745
		9	0.8570
		10	0.8817/,

	BarrelToKT(c),

	Sulphur_GO_nominal(c)
	/	1	0.157
		2	0.293
		3	0.162
		4	0.200
		5	0.263
		6	0.694
		7	0.767
		8	1.550
		9	0.326
		10	1.090/,

	Sulphur_GO_stdev(c),

	VaccuumResidue_nominal(c),

	VaccuumResidue_stdev(c),

	Crude_price(c)
	/	1	115.0
		2	107.5
		3	109.7
		4	110.7
		5	108.4
		6	101.6
		7	114.3
		8	101.3
		9	109.4
		10	104.09/,

	Demand_quantity(p)
	/	1	5
		2	0
		3	0
		4	100
		5	100
		6	0
		7	0/,

	Density_PG98_input(PG98_in)
	/	1	0.58
		2	0.665
		3	0.65
		4	0.77
		5	0.80
		6	0.75/,

	Density_products(p)
	/	1	0.79
		2	0.76
		3	0.75
		4	0.87
		5	0.98
		6	0.54
		7	0.65/,

	Product_VP(p)
	/	1	0.65
		2	0.65
		3	0
		4	0
		5	0
		6	0
		7	0/,

	Product_RON(p)
	/	1	105
		2	100
		3	0
		4	0
		5	0
		6	0
		7	0/,

	Product_MON(p)
	/	1	100
		2	96
		3	0
		4	0
		5	0
		6	0
		7	0/,

	Product_Sulphur(p)
	/	1	0
		2	0
		3	0.001
		4	0
		5	0
		6	0
		7	0/,

	Import_upper(p)
	/	1	0
		2	0
		3	0
		4	0
		5	0
		6	0
		7	0/,

	RON(PG98_in)
	/	1	0
		2	91
		3	71
		4	95
		5	100
		6	93/,

	MON(PG98_in)
	/	1	0
		2	86
		3	68
		4	86
		5	91
		6	82/,

	VP(PG98_in)
	/	1	0
		2	0.4
		3	0.8
		4	0.5
		5	0.5
		6	0.65/,

	HFO_density(c)
	/	1	0.9385
		2	0.9682
		3	0.9423
		4	0.9433
		5	0.9652
		6	0.9727
		7	0.9470
		8	0.9799
		9	0.9562
		10	0.9685/,

	GO_density(c)
	/	1	0.8506
		2	0.8590
		3	0.8413
		4	0.8450
		5	0.8573
		6	0.8688
		7	0.8404
		8	0.8467
		9	0.8477
		10	0.8558/,

	Viscosity_HF1(c)
	/	1	32.5
		2	69.6
		3	38.2
		4	42.7
		5	86.5
		6	75.5
		7	42.3
		8	45.0
		9	53.5
		10	55.2/,

	Viscosity_HF3(c)
	/	1	2.52
		2	2.92
		3	2.61
		4	2.56
		5	2.65
		6	2.95
		7	2.50
		8	2.51
		9	2.62
		10	2.67/,

	Viscosity_products(p)
	/	4	31.5/,

	Sulphur_3(AGO_in)
	/	1	0.1/,

	Crude_lower_bound(c),

	Crude_upper_bound(c),

	Sulphur_2(c,h),

	Sulphur_GO_data(c,h),

	VaccuumResidue_data(c,h),

	Crude_yield_data(c,w,h),

	Desulphurisation_cost(c,h),

	prob(h) "probability of each scenario"
;

loop(c,
	BarrelToKT(c) = (GranularityOfBarrels/6.29)*(Crude_density(c)/1000);
	Sulphur_GO_stdev(c) = 0.1*Sulphur_GO_nominal(c);
	VaccuumResidue_nominal(c) = sum(w$(ord(w)=card(w)), Crude_yield(c,w));
	VaccuumResidue_stdev(c) = 0.1*VaccuumResidue_nominal(c);
);



*=========== Generate scenarios for the uncertain parameters =============

$INCLUDE refinery1_data/5.gms


*======================================================================



loop(h,
	loop(c,
		loop(w$(ord(w) < card(w)),
			Crude_yield_data(c,w,h) = sum(w2$(ord(w2)=card(w2)), Crude_yield(c,w)/(1-Crude_yield(c,w2))*(1-VaccuumResidue_data(c,h)));
		);
		Crude_yield_data(c,w2,h)$(ord(w2)=card(w2)) = VaccuumResidue_data(c,h);
	);
);

loop(c,
	loop(h,
		Desulphurisation_cost(c,h) = ((Sulphur_GO_data(c,h)*109.0909 + 365.4546)/1000)*(0.85/0.159)/GO_density(c);
		Sulphur_2(c,h) = Sulphur_GO_data(c,h)*0.005;
	);
);

loop(c,
	Crude_lower_bound(c) = (Barrel_lower_bound/GranularityOfBarrels)*BarrelToKT(c);
	Crude_upper_bound(c) = (Barrel_upper_bound/GranularityOfBarrels)*BarrelToKT(c);
);





*---------------------------------------------------------------
*-------------    FULL SPACE PROBLEM DEFINITION   --------------
*---------------------------------------------------------------

BINARY VARIABLES
	pickCrude(c)
;

POSITIVE VARIABLES
	crudeQuantity(c)
	flow_Reformer95(h)
	flow_Reformer100(h)
	flow_Cracker_Mogas(h)
	flow_Cracker_AGO(h)
	flow_Isomerisation(h)
	flow_Desulphurisation_CGO(h)
	flow_LG_producing(h)
	flow_LN_producing(h)
	flow_HF_2(h)
	volume_PG98(h)
	volume_ES95(h)
	volume_HF(h)

	blin_CDU_LG(LG_out,h)
	blin_Reformer95_LG(LG_out,h)
	blin_Reformer100_LG(LG_out,h)
	blin_Mogas_LG(LG_out,h)
	blin_AGO_LG(LG_out,h)
	blin_Cracker_Mogas(Cr_CGO,h)
	blin_Cracker_AGO(Cr_CGO,h)

	flow_Desulphurisation_1(c,h)
	flow_AGO_1(c,h)
	flow_AGO_2(c,h)
	flow_HF_1(c,h)
	flow_HF_3(c,h)
	flow_Burn(Burn,h)
	flow_PG98(PG98_in,h)
	flow_ES95(PG98_in,h)
	flow_AGO_3(AGO_in,h)
	flow_JPF(JPF_out,h)
	flow_Import(p,h)
	fraction_LG(LG_in,h)
	fraction_CGO(Cr_mode,h)
;

VARIABLES
	objvar
;

EQUATIONS
	Mass_balance1(h)
	Mass_balance2(h)
	Mass_balance3(h)
	Mass_balance4(h)
	Mass_balance5(h)
	Mass_balance7(h)
	GO_balance(c,h)
	VR_balance(c,h)
	Desulphurisation_balance(c,h)
	Demand_constraint1(h)
	Demand_constraint2(h)
	Demand_constraint3(h)
	Demand_constraint4(h)
	Demand_constraint5(h)
	Demand_constraint6(h)
	Demand_constraint7(h)

	CDU_capacity_bound
	Crude_bound(c)
	Crude_selection(c)
	Desulphurisation_capacity_bound(h)
	Cracker_capacity_bound(h)
	Reformer_capacity_bound(h)
	Reformer95_balance(h)
	Reformer100_Balance(h)
	Isomerisation_balance(h)
	CN_balance(h)
	CGO_balance(h)
	Desulphurisation_CGO_balance(h)
	PG98_volume_def(h)
	ES95_volume_def(h)
	Butane95_constraint(h)
	Butane98_constraint(h)
	LG_split_balance(h)
	LG_balance(h)
	Reformer95_LG_balance(h)
	Reformer100_LG_balance(h)
	Cracker_Mogas_LG_balance(h)
	Cracker_AGO_LG_balance(h)
	pq_ES95_constraint(h)
	pq_PG98_constraint(h)
	pq_demand_constraint(h)
	pq_burn_constraint(h)
	VP_ES95_upper(h)
	VP_ES95_lower(h)
	VP_PG98_upper(h)
	VP_PG98_lower(h)
	RON_PG98(h)
	RON_ES95(h)
	Sensitivity_PG98(h)
	Sensitivity_ES95(h)

	blincon_CDU_LG1(h)
	blincon_CDU_LG2(h)
	blincon_CDU_LG3(h)
	blincon_CDU_LG4(h)

	blincon_Reformer95_LG1(h)
	blincon_Reformer95_LG2(h)
	blincon_Reformer95_LG3(h)
	blincon_Reformer95_LG4(h)

	blincon_Reformer100_LG1(h)
	blincon_Reformer100_LG2(h)
	blincon_Reformer100_LG3(h)
	blincon_Reformer100_LG4(h)

	blincon_Mogas_LG1(h)
	blincon_Mogas_LG2(h)
	blincon_Mogas_LG3(h)
	blincon_Mogas_LG4(h)

	blincon_AGO_LG1(h)
	blincon_AGO_LG2(h)
	blincon_AGO_LG3(h)
	blincon_AGO_LG4(h)

	blincon_Cracker_Mogas1(h)
	blincon_Cracker_Mogas2(h)
	blincon_Cracker_Mogas3(h)

	blincon_Cracker_AGO1(h)
	blincon_Cracker_AGO2(h)
	blincon_Cracker_AGO3(h)

	pq_AGO_constraint(h)
	pq_HF_constraint(h)
	pq_Desulphurisation_constraint(h)
	CGO_split_balance(h)
	Cracker_Mogas_CGO_balance(h)
	Cracker_AGO_CGO_balance(h)
	HF_volume_def(h)
	AGO_sulphur_balance(h)
	HF_viscosity_upper(h)
	HF_viscosity_lower(h)
	Refinery_Fuel(h)

 	objfn
;

*-------------------------------------
*		EQUATION DEFINITIONS
*-------------------------------------

	CDU_capacity_bound .. 	sum(c, crudeQuantity(c)*BarrelToKT(c)/GranularityOfBarrels) =l= CDU_capacity
;

	Crude_selection(c) .. crudeQuantity(c) =g= pickCrude(c)*Barrel_lower_bound
;

	Crude_bound(c) .. 	CrudeQuantity(c) =l= pickCrude(c)*Barrel_upper_bound
;

	Desulphurisation_capacity_bound(h) .. flow_Desulphurisation_CGO(h) + sum(c, flow_Desulphurisation_1(c,h)) =l= Desulphurisation_capacity
;

	Mass_balance1(h) .. Reformer_fraction('1','1')*flow_Reformer95(h) +
						Reformer_fraction('2','1')*flow_Reformer100(h) +
						Cracker_fraction('1','1')*flow_Cracker_Mogas(h) +
						Cracker_fraction('2','1')*flow_Cracker_AGO(h) +
						Isomerisation_fraction('1')*flow_Isomerisation(h) +
						Desulphurisation_fraction2('2')*flow_Desulphurisation_CGO(h) -
						flow_Burn('1',h) +
						sum(c,
							Crude_yield_data(c,'1',h)*crudeQuantity(c)*BarrelToKT(c)/GranularityOfBarrels +
							Desulphurisation_fraction(c,'2')*flow_Desulphurisation_1(c,h)
						) =e= 0
;

	Mass_balance2(h) .. Reformer_fraction('1','2')*flow_Reformer95(h) +
						Reformer_fraction('2','2')*flow_Reformer100(h) +
						Cracker_fraction('1','2')*flow_Cracker_Mogas(h) +
						Cracker_fraction('2','2')*flow_Cracker_AGO(h) -
						flow_LG_producing(h) - flow_PG98('1',h) -
						flow_ES95('1',h) - flow_Burn('2',h) +
						sum(c,
							Crude_yield_data(c,'2',h)*crudeQuantity(c)*BarrelToKT(c)/GranularityOfBarrels
						) =e= 0
;

	Mass_balance3(h) .. -flow_LN_producing(h) - flow_Burn('3',h) -
						flow_PG98('3',h) - flow_ES95('3',h) -
						flow_Isomerisation(h) - flow_JPF('1',h)*JPF_fraction('1','1') -
						flow_JPF('2',h)*JPF_fraction('1','2') +
						sum(c,
							Crude_yield_data(c,'3',h)*crudeQuantity(c)*BarrelToKT(c)/GranularityOfBarrels
						) =e= 0
;

	Mass_balance4(h) .. -flow_JPF('1',h)*JPF_fraction('2','1') -
						flow_JPF('2',h)*JPF_fraction('2','2') -
						flow_Reformer95(h) - flow_Reformer100(h) +
						sum(c,
							Crude_yield_data(c,'4',h)*crudeQuantity(c)*BarrelToKT(c)/GranularityOfBarrels
						) =e= 0
;

	Mass_balance5(h) .. -flow_JPF('1',h)*JPF_fraction('3','1') -
						flow_JPF('2',h)*JPF_fraction('3','2') -
						flow_AGO_3('1',h) +
						sum(c,
							Crude_yield_data(c,'5',h)*crudeQuantity(c)*BarrelToKT(c)/GranularityOfBarrels
						) =e= 0
;

	Mass_balance7(h) .. -flow_Cracker_Mogas(h) - flow_Cracker_AGO(h) +
						sum(c,
							Crude_yield_data(c,'7',h)*crudeQuantity(c)*BarrelToKT(c)/GranularityOfBarrels
						) =e= 0
;

	GO_balance(c,h) .. -flow_AGO_1(c,h) - flow_Desulphurisation_1(c,h) - flow_HF_3(c,h) +
						Crude_yield_data(c,'6',h)*crudeQuantity(c)*BarrelToKT(c)/GranularityOfBarrels =e= 0
;

	VR_balance(c,h) .. 	Crude_yield_data(c,'8',h)*crudeQuantity(c)*BarrelToKT(c)/GranularityOfBarrels =e= flow_HF_1(c,h)
;

	Desulphurisation_balance(c,h) .. Desulphurisation_fraction(c,'1')*flow_Desulphurisation_1(c,h) =e= flow_AGO_2(c,h)
;

	Reformer95_balance(h) .. 	flow_Reformer95(h)*Reformer_fraction('1','3') +
								flow_Reformer100(h)*Reformer_fraction('2','3') =e=
								flow_PG98('4',h) + flow_ES95('4',h)
;

	Reformer100_balance(h) .. 	flow_Reformer95(h)*Reformer_fraction('1','4') +
								flow_Reformer100(h)*Reformer_fraction('2','4') =e=
								flow_PG98('5',h) + flow_ES95('5',h)
;

	Isomerisation_balance(h) .. flow_Isomerisation(h)*Isomerisation_fraction('2') =e=
								flow_PG98('2',h) + flow_ES95('2',h)
;

	CN_balance(h) .. 	flow_Cracker_Mogas(h)*Cracker_fraction('1','3') +
						flow_Cracker_AGO(h)*Cracker_fraction('2','3') =e=
						flow_PG98('6',h) + flow_ES95('6',h)
;

	CGO_balance(h) .. 	flow_Cracker_Mogas(h)*Cracker_fraction('1','4') +
						flow_Cracker_AGO(h)*Cracker_fraction('2','4') =e=
						flow_Desulphurisation_CGO(h) + flow_HF_2(h) + flow_AGO_3('2',h)
;

	Desulphurisation_CGO_balance(h) .. Desulphurisation_fraction2('1')*flow_Desulphurisation_CGO(h) =e= flow_AGO_3('3',h)
;

	Demand_constraint1(h) .. flow_Import('1',h) + sum(PG98_in, flow_PG98(PG98_in,h)) =g= Demand_quantity('1')
;

	Demand_constraint2(h) .. flow_Import('2',h) + sum(PG98_in, flow_ES95(PG98_in,h)) =g= Demand_quantity('2')
;

	Demand_constraint3(h) .. flow_Import('3',h) + sum(JPF_out, flow_JPF(JPF_out,h)) =g= Demand_quantity('3')
;

	Demand_constraint4(h) .. flow_Import('4',h) + sum(AGO_in, flow_AGO_3(AGO_in,h)) +
								sum(c, flow_AGO_1(c,h) + flow_AGO_2(c,h)) =g= Demand_quantity('4')
;

	Demand_constraint5(h) .. flow_Import('5',h) + flow_HF_2(h) + 
								sum(c, flow_HF_1(c,h) + flow_HF_3(c,h)) =g= Demand_quantity('5')
;

	Demand_constraint6(h) .. flow_Import('6',h) + flow_LG_producing(h) =g= Demand_quantity('6')
;

	Demand_constraint7(h) .. flow_Import('7',h) + flow_LN_producing(h) =g= Demand_quantity('7')
;

	PG98_volume_def(h) ..	flow_Import('1',h)/Density_products('1') +
							sum(PG98_in, flow_PG98(PG98_in,h)/Density_PG98_input(PG98_in)) =e= volume_PG98(h)
;

	ES95_volume_def(h) ..	flow_Import('2',h)/Density_products('2') +
							sum(PG98_in, flow_ES95(PG98_in,h)/Density_PG98_input(PG98_in)) =e= volume_ES95(h)
;

	Butane95_constraint(h) ..	flow_ES95('1',h)/Density_PG98_input('1') +
								0.03*flow_Import('2',h)/Density_products('2') =l= 0.05*volume_ES95(h)
;

	Butane98_constraint(h) ..	flow_PG98('1',h)/Density_PG98_input('1') +
								0.03*flow_Import('1',h)/Density_products('2') =l= 0.05*volume_PG98(h)
;

	blincon_CDU_LG1(h) .. blin_CDU_LG('1',h) =e= fraction_LG('1',h)*flow_ES95('1',h)
;

	blincon_CDU_LG2(h) .. blin_CDU_LG('2',h) =e= fraction_LG('1',h)*flow_PG98('1',h)
;

	blincon_CDU_LG3(h) .. blin_CDU_LG('3',h) =e= fraction_LG('1',h)*flow_Burn('2',h)
;

	blincon_CDU_LG4(h) .. blin_CDU_LG('4',h) =e= fraction_LG('1',h)*flow_LG_producing(h)
;

	blincon_Reformer95_LG1(h) .. blin_Reformer95_LG('1',h) =e= fraction_LG('2',h)*flow_ES95('1',h)
;

	blincon_Reformer95_LG2(h) .. blin_Reformer95_LG('2',h) =e= fraction_LG('2',h)*flow_PG98('1',h)
;

	blincon_Reformer95_LG3(h) .. blin_Reformer95_LG('3',h) =e= fraction_LG('2',h)*flow_Burn('2',h)
;

	blincon_Reformer95_LG4(h) .. blin_Reformer95_LG('4',h) =e= fraction_LG('2',h)*flow_LG_producing(h)
;

	blincon_Reformer100_LG1(h) .. blin_Reformer100_LG('1',h) =e= fraction_LG('3',h)*flow_ES95('1',h)
;

	blincon_Reformer100_LG2(h) .. blin_Reformer100_LG('2',h) =e= fraction_LG('3',h)*flow_PG98('1',h)
;

	blincon_Reformer100_LG3(h) .. blin_Reformer100_LG('3',h) =e= fraction_LG('3',h)*flow_Burn('2',h)
;

	blincon_Reformer100_LG4(h) .. blin_Reformer100_LG('4',h) =e= fraction_LG('3',h)*flow_LG_producing(h)
;

	blincon_Mogas_LG1(h) .. blin_Mogas_LG('1',h) =e= fraction_LG('4',h)*flow_ES95('1',h)
;

	blincon_Mogas_LG2(h) .. blin_Mogas_LG('2',h) =e= fraction_LG('4',h)*flow_PG98('1',h)
;

	blincon_Mogas_LG3(h) .. blin_Mogas_LG('3',h) =e= fraction_LG('4',h)*flow_Burn('2',h)
;

	blincon_Mogas_LG4(h) .. blin_Mogas_LG('4',h) =e= fraction_LG('4',h)*flow_LG_producing(h)
;

	blincon_AGO_LG1(h) .. blin_AGO_LG('1',h) =e= fraction_LG('5',h)*flow_ES95('1',h)
;

	blincon_AGO_LG2(h) .. blin_AGO_LG('2',h) =e= fraction_LG('5',h)*flow_PG98('1',h)
;

	blincon_AGO_LG3(h) .. blin_AGO_LG('3',h) =e= fraction_LG('5',h)*flow_Burn('2',h)
;

	blincon_AGO_LG4(h) .. blin_AGO_LG('4',h) =e= fraction_LG('5',h)*flow_LG_producing(h)
;

	LG_balance(h) .. 	sum(LG_out, blin_CDU_LG(LG_out,h)) =e=
						sum(c,
							Crude_yield_data(c,'2',h)*crudeQuantity(c)*BarrelToKT(c)/GranularityOfBarrels
						)
;

	Reformer95_LG_balance(h) .. flow_Reformer95(h)*Reformer_fraction('1','2') =e=
								sum(LG_out, blin_Reformer95_LG(LG_out,h))
;

	Reformer100_LG_balance(h) .. 	flow_Reformer100(h)*Reformer_fraction('2','2') =e=
									sum(LG_out, blin_Reformer100_LG(LG_out,h))
;

	Cracker_Mogas_LG_balance(h) .. 	flow_Cracker_Mogas(h)*Cracker_fraction('1','2') =e=
									sum(LG_out, blin_Mogas_LG(LG_out,h))
;

	Cracker_AGO_LG_balance(h) .. 	flow_Cracker_AGO(h)*Cracker_fraction('2','2') =e=
									sum(LG_out, blin_AGO_LG(LG_out,h))
;

	pq_ES95_constraint(h) .. 	blin_CDU_LG('1',h) + blin_Reformer95_LG('1',h) +
								blin_Reformer100_LG('1',h) + blin_Mogas_LG('1',h) +
								blin_AGO_LG('1',h) =e= flow_ES95('1',h)
;

	pq_PG98_constraint(h) .. 	blin_CDU_LG('2',h) + blin_Reformer95_LG('2',h) +
								blin_Reformer100_LG('2',h) + blin_Mogas_LG('2',h) +
								blin_AGO_LG('2',h) =e= flow_PG98('1',h)
;

	pq_burn_constraint(h) .. 	blin_CDU_LG('3',h) + blin_Reformer95_LG('3',h) +
								blin_Reformer100_LG('3',h) + blin_Mogas_LG('3',h) +
								blin_AGO_LG('3',h) =e= flow_Burn('2',h)
;

	pq_demand_constraint(h) .. 	blin_CDU_LG('4',h) + blin_Reformer95_LG('4',h) +
								blin_Reformer100_LG('4',h) + blin_Mogas_LG('4',h) +
								blin_AGO_LG('4',h) =e= flow_LG_producing(h)
;

	LG_split_balance(h) .. sum(LG_in, fraction_LG(LG_in,h)) =e= 1
;

	VP_ES95_lower(h) .. -0.45*volume_ES95(h) + flow_Import('2',h)*Product_VP('2')/Density_products('2') +
						sum(PG98_in, VP(PG98_in)*flow_ES95(PG98_in,h)/Density_PG98_input(PG98_in)) +
						LG_parameters('1','1')*blin_CDU_LG('1',h)/Density_PG98_input('1') +
						LG_parameters('1','2')*blin_Reformer95_LG('1',h)/Density_PG98_input('1') +
						LG_parameters('1','3')*blin_Reformer100_LG('1',h)/Density_PG98_input('1') +
						LG_parameters('1','4')*blin_Mogas_LG('1',h)/Density_PG98_input('1') +
						LG_parameters('1','5')*blin_AGO_LG('1',h)/Density_PG98_input('1') =g= 0
;

	VP_ES95_upper(h) .. -0.80*volume_ES95(h) + flow_Import('2',h)*Product_VP('2')/Density_products('2') +
						sum(PG98_in, VP(PG98_in)*flow_ES95(PG98_in,h)/Density_PG98_input(PG98_in)) +
						LG_parameters('1','1')*blin_CDU_LG('1',h)/Density_PG98_input('1') +
						LG_parameters('1','2')*blin_Reformer95_LG('1',h)/Density_PG98_input('1') +
						LG_parameters('1','3')*blin_Reformer100_LG('1',h)/Density_PG98_input('1') +
						LG_parameters('1','4')*blin_Mogas_LG('1',h)/Density_PG98_input('1') +
						LG_parameters('1','5')*blin_AGO_LG('1',h)/Density_PG98_input('1') =l= 0
;

	VP_PG98_lower(h) .. -0.50*volume_PG98(h) + flow_Import('1',h)*Product_VP('1')/Density_products('1') +
						sum(PG98_in, VP(PG98_in)*flow_PG98(PG98_in,h)/Density_PG98_input(PG98_in)) +
						LG_parameters('1','1')*blin_CDU_LG('2',h)/Density_PG98_input('1') +
						LG_parameters('1','2')*blin_Reformer95_LG('2',h)/Density_PG98_input('1') +
						LG_parameters('1','3')*blin_Reformer100_LG('2',h)/Density_PG98_input('1') +
						LG_parameters('1','4')*blin_Mogas_LG('2',h)/Density_PG98_input('1') +
						LG_parameters('1','5')*blin_AGO_LG('2',h)/Density_PG98_input('1') =g= 0
;

	VP_PG98_upper(h) .. -0.86*volume_PG98(h) + flow_Import('1',h)*Product_VP('1')/Density_products('1') +
						sum(PG98_in, VP(PG98_in)*flow_PG98(PG98_in,h)/Density_PG98_input(PG98_in)) +
						LG_parameters('1','1')*blin_CDU_LG('2',h)/Density_PG98_input('1') +
						LG_parameters('1','2')*blin_Reformer95_LG('2',h)/Density_PG98_input('1') +
						LG_parameters('1','3')*blin_Reformer100_LG('2',h)/Density_PG98_input('1') +
						LG_parameters('1','4')*blin_Mogas_LG('2',h)/Density_PG98_input('1') +
						LG_parameters('1','5')*blin_AGO_LG('2',h)/Density_PG98_input('1') =l= 0
;

	RON_PG98(h) .. 	-98*volume_PG98(h) + flow_Import('1',h)*Product_RON('1')/Density_products('1') +
					sum(PG98_in, RON(PG98_in)*flow_PG98(PG98_in,h)/Density_PG98_input(PG98_in)) +
					LG_parameters('2','1')*blin_CDU_LG('2',h)/Density_PG98_input('1') +
					LG_parameters('2','2')*blin_Reformer95_LG('2',h)/Density_PG98_input('1') +
					LG_parameters('2','3')*blin_Reformer100_LG('2',h)/Density_PG98_input('1') +
					LG_parameters('2','4')*blin_Mogas_LG('2',h)/Density_PG98_input('1') +
					LG_parameters('2','5')*blin_AGO_LG('2',h)/Density_PG98_input('1') =g= 0
;

	RON_ES95(h) .. 	-95*volume_ES95(h) + flow_Import('2',h)*Product_RON('2')/Density_products('2') +
					sum(PG98_in, RON(PG98_in)*flow_ES95(PG98_in,h)/Density_PG98_input(PG98_in)) +
					LG_parameters('2','1')*blin_CDU_LG('1',h)/Density_PG98_input('1') +
					LG_parameters('2','2')*blin_Reformer95_LG('1',h)/Density_PG98_input('1') +
					LG_parameters('2','3')*blin_Reformer100_LG('1',h)/Density_PG98_input('1') +
					LG_parameters('2','4')*blin_Mogas_LG('1',h)/Density_PG98_input('1') +
					LG_parameters('2','5')*blin_AGO_LG('1',h)/Density_PG98_input('1') =g= 0
;

	Sensitivity_PG98(h) .. 	-10*volume_PG98(h) + flow_Import('1',h)*(Product_RON('1') - Product_MON('1'))/Density_products('1') +
							sum(PG98_in, (RON(PG98_in) - MON(PG98_in))*flow_PG98(PG98_in,h)/Density_PG98_input(PG98_in)) +
							(LG_parameters('2','1') - LG_parameters('3','1'))*blin_CDU_LG('2',h)/Density_PG98_input('1') +
							(LG_parameters('2','2') - LG_parameters('3','2'))*blin_Reformer95_LG('2',h)/Density_PG98_input('1') +
							(LG_parameters('2','3') - LG_parameters('3','3'))*blin_Reformer100_LG('2',h)/Density_PG98_input('1') +
							(LG_parameters('2','4') - LG_parameters('3','4'))*blin_Mogas_LG('2',h)/Density_PG98_input('1') +
							(LG_parameters('2','5') - LG_parameters('3','5'))*blin_AGO_LG('2',h)/Density_PG98_input('1') =l= 0
;

	Sensitivity_ES95(h) .. 	-10*volume_ES95(h) + flow_Import('2',h)*(Product_RON('2') - Product_MON('2'))/Density_products('2') +
							sum(PG98_in, (RON(PG98_in) - MON(PG98_in))*flow_ES95(PG98_in,h)/Density_PG98_input(PG98_in)) +
							(LG_parameters('2','1') - LG_parameters('3','1'))*blin_CDU_LG('1',h)/Density_PG98_input('1') +
							(LG_parameters('2','2') - LG_parameters('3','2'))*blin_Reformer95_LG('1',h)/Density_PG98_input('1') +
							(LG_parameters('2','3') - LG_parameters('3','3'))*blin_Reformer100_LG('1',h)/Density_PG98_input('1') +
							(LG_parameters('2','4') - LG_parameters('3','4'))*blin_Mogas_LG('1',h)/Density_PG98_input('1') +
							(LG_parameters('2','5') - LG_parameters('3','5'))*blin_AGO_LG('1',h)/Density_PG98_input('1') =l= 0
;

	blincon_Cracker_Mogas1(h) .. blin_Cracker_Mogas('1',h) =e= fraction_CGO('1',h)*flow_AGO_3('2',h)
;

	blincon_Cracker_Mogas2(h) .. blin_Cracker_Mogas('2',h) =e= fraction_CGO('1',h)*flow_HF_2(h)
;

	blincon_Cracker_Mogas3(h) .. blin_Cracker_Mogas('3',h) =e= fraction_CGO('1',h)*flow_Desulphurisation_CGO(h)
;

	blincon_Cracker_AGO1(h) .. blin_Cracker_AGO('1',h) =e= fraction_CGO('2',h)*flow_AGO_3('2',h)
;

	blincon_Cracker_AGO2(h) .. blin_Cracker_AGO('2',h) =e= fraction_CGO('2',h)*flow_HF_2(h)
;

	blincon_Cracker_AGO3(h) .. blin_Cracker_AGO('3',h) =e= fraction_CGO('2',h)*flow_Desulphurisation_CGO(h)
;

	Cracker_Mogas_CGO_balance(h) .. blin_Cracker_Mogas('1',h) + blin_Cracker_Mogas('2',h) +
									blin_Cracker_Mogas('3',h) =e= flow_Cracker_Mogas(h)*Cracker_fraction('1','4')
;

	Cracker_AGO_CGO_balance(h) .. 	blin_Cracker_AGO('1',h) + blin_Cracker_AGO('2',h) +
									blin_Cracker_AGO('3',h) =e= flow_Cracker_AGO(h)*Cracker_fraction('2','4')
;

	CGO_split_balance(h) .. sum(Cr_mode, fraction_CGO(Cr_mode,h)) =e= 1
;

	pq_AGO_constraint(h) .. blin_Cracker_Mogas('1',h) + blin_Cracker_AGO('1',h) =e= flow_AGO_3('2',h)
;

	pq_HF_constraint(h) .. blin_Cracker_Mogas('2',h) + blin_Cracker_AGO('2',h) =e= flow_HF_2(h)
;

	pq_Desulphurisation_constraint(h) .. blin_Cracker_Mogas('3',h) + blin_Cracker_AGO('3',h) =e= flow_Desulphurisation_CGO(h)
;

	HF_volume_def(h) .. -volume_HF(h) + flow_Import('5',h)/Density_products('5') +
						flow_HF_2(h)/CGO_density +
						sum(c, flow_HF_1(c,h)/HFO_density(c) + flow_HF_3(c,h)/GO_density(c)) =e= 0
;

	HF_viscosity_lower(h) .. 	flow_Import('5',h)*Viscosity_products('5')/Density_products('5') +
								sum(c,
									flow_HF_1(c,h)*Viscosity_HF1(c)/HFO_density(c) +
									flow_HF_3(c,h)*Viscosity_HF3(c)/GO_density(c)
								) +
								(blin_Cracker_Mogas('2',h)*Mogas_viscosity + blin_Cracker_AGO('2',h)*AGO_viscosity)/CGO_density -
								30*volume_HF(h) =g= 0
;

	HF_viscosity_upper(h) .. 	flow_Import('5',h)*Viscosity_products('5')/Density_products('5') +
								sum(c,
									flow_HF_1(c,h)*Viscosity_HF1(c)/HFO_density(c) +
									flow_HF_3(c,h)*Viscosity_HF3(c)/GO_density(c)
								) +
								(blin_Cracker_Mogas('2',h)*Mogas_viscosity + blin_Cracker_AGO('2',h)*AGO_viscosity)/CGO_density -
								33*volume_HF(h) =l= 0
;

	AGO_sulphur_balance(h) .. 	flow_Import('4',h)*Product_sulphur('4') - Sulphur_spec*flow_Import('4',h) +
								sum(c,
									(Sulphur_GO_data(c,h) - Sulphur_spec)*flow_AGO_1(c,h) +
									(Sulphur_2(c,h) - Sulphur_spec)*flow_AGO_2(c,h)
								) +
								flow_AGO_3('1',h)*(Sulphur_3('1') - Sulphur_spec) +
								blin_Cracker_AGO('1',h)*(AGO_sulphur - Sulphur_spec) +
								blin_Cracker_Mogas('1',h)*(Mogas_sulphur - Sulphur_spec) +
								blin_Cracker_AGO('3',h)*AGO_sulphur*0.005 +
								blin_Cracker_Mogas('3',h)*Mogas_sulphur*0.005 -
								Sulphur_spec*flow_AGO_3('3',h) =l= 0
;

	Refinery_Fuel(h) .. 1.3*flow_Burn('1',h) + 1.2*flow_Burn('2',h) + 1.1*flow_Burn('3',h) -
						flow_Reformer95(h)*Reformer_fraction('1','5') -
						flow_Reformer100(h)*Reformer_fraction('2','5') -
						flow_Cracker_Mogas(h)*Cracker_fraction('1','5') -
						flow_Cracker_AGO(h)*Cracker_fraction('2','5') -
						flow_Isomerisation(h)*Isomerisation_fraction('3') -
						flow_Desulphurisation_CGO(h)*Desulphurisation_fraction2('3') - 15.2 -
						sum(c,
							0.018*crudeQuantity(c)*BarrelToKT(c)/GranularityOfBarrels +
							flow_Desulphurisation_1(c,h)*Desulphurisation_fraction(c,'3')
						) =g= 0
;

	Cracker_capacity_bound(h) .. flow_Cracker_Mogas(h) + flow_Cracker_AGO(h) =l= Cracker_capacity
;

	Reformer_capacity_bound(h) .. flow_Reformer95(h) + flow_Reformer100(h) =l= Reformer_capacity
;


    objfn .. objvar =e= sum(h,prob(h)*(
								Cracker_Mogas_cost*flow_Cracker_Mogas(h) +
								Cracker_AGO_cost*flow_Cracker_AGO(h) +
								Reformer95_cost*flow_Reformer95(h) +
								Reformer100_cost*flow_Reformer100(h) +
								Isomerisation_cost*flow_Isomerisation(h) +
								Desulphurisation_CGO_cost*flow_Desulphurisation_CGO(h) -
								LG_sale*flow_LG_producing(h) -
								LN_sale*flow_LN_producing(h) -
								HF_sale*flow_HF_2(h) +
								sum(c,
									Desulphurisation_cost(c,h)*flow_Desulphurisation_1(c,h) -
									AGO_sale*flow_AGO_1(c,h) -
									AGO_sale*flow_AGO_2(c,h) -
									HF_sale*flow_HF_1(c,h) -
									HF_sale*flow_HF_3(c,h) +
									(crudeQuantity(c)/1000)*(Crude_price(c)+1)
								) -
								sum(PG98_in,
									PG98_sale*flow_PG98(PG98_in,h) +
									ES95_sale*flow_ES95(PG98_in,h)
								) -
								sum(JPF_out,
									JET_sale*flow_JPF(JPF_out,h)
								) -
								sum(AGO_in,
									AGO_sale*flow_AGO_3(AGO_in,h)
								)
							) );


MODEL refinery1 /all/;


*---------------------------------------------
*			SET VARIABLE BOUNDS
*---------------------------------------------

flow_Desulphurisation_1.up(c,h) = Desulphurisation_capacity;
flow_AGO_1.up(c,h) = Crude_upper_bound(c);
flow_AGO_2.up(c,h) = Desulphurisation_capacity;
flow_HF_1.up(c,h) = Crude_upper_bound(c);
flow_HF_3.up(c,h) = Crude_upper_bound(c);
flow_PG98.up(PG98_in,h) = CDU_capacity;
flow_ES95.up(PG98_in,h) = CDU_capacity;
flow_Burn.up(Burn,h) = CDU_capacity;
flow_AGO_3.up(AGO_in,h) = Cracker_capacity;
flow_AGO_3.up('1',h) = CDU_capacity;
flow_JPF.up(JPF_out,h) = CDU_capacity;
flow_Import.up(p,h) = Import_upper(p);
fraction_LG.up(LG_in,h) = 1;
flow_Reformer95.lo(h) = Reformer95_lower;
flow_Reformer95.up(h) = Reformer_capacity;
flow_Reformer100.up(h) = Reformer_capacity - Reformer95_lower;
flow_Cracker_Mogas.up(h) = Cracker_capacity;
flow_Cracker_AGO.up(h) = Cracker_capacity;
flow_Isomerisation.up(h) = CDU_capacity;
flow_Desulphurisation_CGO.up(h) = Cracker_capacity;
flow_LG_producing.up(h) = CDU_capacity;
flow_LN_producing.up(h) = CDU_capacity;
flow_HF_2.up(h) = Cracker_capacity;
volume_PG98.up(h) = CDU_capacity/Density_PG98_input('1');
volume_ES95.up(h) = CDU_capacity/Density_PG98_input('1');
fraction_CGO.up(Cr_mode,h) = 1;
volume_HF.up(h) = CDU_capacity/GO_density('7');
blin_CDU_LG.up(LG_out,h) = CDU_capacity;
blin_Reformer95_LG.up(LG_out,h) = CDU_capacity;
blin_Reformer100_LG.up(LG_out,h) = CDU_capacity;
blin_Mogas_LG.up(LG_out,h) = CDU_capacity;
blin_AGO_LG.up(LG_out,h) = CDU_capacity;
blin_Cracker_Mogas.up(Cr_CGO,h) = Cracker_capacity;
blin_Cracker_AGO.up(Cr_CGO,h) = Cracker_Capacity;



*------------------------------------------
*			SOLVE THE PROBLEM
*------------------------------------------

solve refinery1 minimizing objvar using MINLP;


*---------------------------------------------------
*			PRINT FINAL SOLUTION STATISTICS
*---------------------------------------------------

SCALARS solver_time, wall_time;

solver_time = refinery1.resusd;
wall_time = refinery1.etsolve;

display solver_time, wall_time;

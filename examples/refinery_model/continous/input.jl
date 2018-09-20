using JuMP
using BARON 
using CPLEX

#define sets 
crudes = 1:10 
components = 1:8 
Re_in= 1:2 
Re_out= 1:6 
Cr_in= 1:2 
Cr_out= 1:6 
Cr_CGO= 1:3 
Cr_mode= 1:2 
Iso_out= 1:4 
Des_out= 1:4 
PG98_in= 1:6 
Burn= 1:3 
JPF_in= 1:3 
JPF_out= 1:2 
AGO_in= 1:3 
products = 1:7 
LG_in= 1:5 
LG_out= 1:4 
LG_prop= 1:3 

#scalars
Desulphurisation_capacity  = 125 
CDU_capacity  = 700 
Reformer95_lower  = 5 
Reformer_capacity  = 65 
Cracker_capacity  = 175 
GranularityOfBarrels  = 5000 
LG_sale  = 561.6 
LN_sale  = 1003 
HF_sale  = 637 
ES95_sale  = 1194 
PG98_sale  = 1231 
JET_sale  = 923 
AGO_sale  = 907 
CGO_density  = 0.95 
Mogas_viscosity  = 12.2 
AGO_viscosity  = 11.65 
Mogas_Sulphur  = 2.1 
AGO_Sulphur  = 1.68 

Isomerisation_cost  = 6 
Reformer95_cost  = 2.7 
Reformer100_cost  = 3.2 
Cracker_Mogas_cost  = 3.2 
Cracker_AGO_cost  = 3 
Barrel_lower_bound  = 100000 
Barrel_upper_bound  = 1500000 
Sulphur_spec  = 0.0015 

Desulphurisation_CGO_cost = ((Mogas_Sulphur*109.0909 + 365.4546)/1000)*(0.85/0.159)/CGO_density


#tables

Reformer_fraction =[	
	0.08  0.09  0.83  0     0.019 2.7
	0.09  0.12  0     0.79  0.026 3.2]


Cracker_fraction =[	
	0.015 0.053 0.436 0.446 0.007 3.2
	0.012 0.046 0.381 0.511 0.007 3.0]


Desulphurisation_fraction =[	
	0.98  0.02  0.02  0
	0.98  0.02  0.02  0
	0.98  0.02  0.02  0
	0.98  0.02  0.02  0
	0.98  0.02  0.02  0
	0.97  0.03  0.02  0
	0.97  0.03  0.02  0
	0.96  0.04  0.02  0
	0.98  0.02  0.02  0
		0.96  0.04  0.02  0]


JPF_fraction =[	
	0.05  0.035
	0.10  0.065
	0.85  0.900]


Crude_yield =[	
	0.0020 0.0091 0.0698 0.1598 0.1003 0.2876 0.2682 0.1032
	0.0020 0.0089 0.0480 0.0959 0.0796 0.2249 0.2735 0.2672
	0.0020 0.0080 0.0610 0.1206 0.0861 0.2414 0.2646 0.2163
	0.0040 0.0200 0.0851 0.1532 0.0947 0.2539 0.2535 0.1356
	0.0020 0.0115 0.0543 0.1026 0.0765 0.2286 0.2695 0.2550
	0.0010 0.0064 0.0246 0.0607 0.0518 0.1900 0.2932 0.3723
	0.0020 0.0155 0.0945 0.1661 0.1160 0.2656 0.2317 0.1086
	0.0029 0.0130 0.0652 0.1196 0.0838 0.2127 0.2408 0.2620
	0.0040 0.0157 0.0749 0.1267 0.0915 0.2353 0.2510 0.2009
		0.0040 0.0107 0.0604 0.1123 0.0784 0.2092 0.2491 0.2759]


LG_parameters =[	
	4.30  4.28  4.36  4.21  4.22
	93.0  92.5  93.6  92.7  93.9
	90.0  89.6  90.9  89.0  90.2]

#parameters 

Isomerisation_fraction=[	0.03
0.97
0.04
6.0]

Desulphurisation_fraction2=[	0.96
0.04
0.02
20.0]

Crude_density=[	0.8441
0.8910
0.8441
0.8369
0.8829
0.9315
0.8252
0.8745
0.8570
	0.8817]



Sulphur_GO_nominal=[	0.157
0.293
0.162
0.200
0.263
0.694
0.767
1.550
0.326
	1.090]

Crude_price=[	115.0
107.5
109.7
110.7
108.4
101.6
114.3
101.3
109.4
	104.09]

Demand_quantity=[	5
0
0
100
100
0
0]

Density_PG98_input=[	0.58
0.665
0.65
0.77
0.80
0.75]

Density_products=[	0.79
0.76
0.75
0.87
0.98
0.54
0.65]

Product_VP=[	0.65
0.65
0
0
0
0
0]

Product_RON=[	105
100
0
0
0
0
0]

Product_MON=[	100
96
0
0
0
0
0]

Product_Sulphur=[	0
0
0.001
0
0
0
0]

Import_upper=[	0
0
0
0
0
0
0]

RON=[	0
91
71
95
100
93]

MON=[	0
86
68
86
91
82]

VP=[	0
0.4
0.8
0.5
0.5
0.65]

HFO_density=[	0.9385
0.9682
0.9423
0.9433
0.9652
0.9727
0.9470
0.9799
0.9562
	0.9685]

GO_density=[	0.8506
0.8590
0.8413
0.8450
0.8573
0.8688
0.8404
0.8467
0.8477
	0.8558]

Viscosity_HF1=[	32.5
69.6
38.2
42.7
86.5
75.5
42.3
45.0
53.5
	55.2]

Viscosity_HF3=[	2.52
2.92
2.61
2.56
2.65
2.95
2.50
2.51
2.62
	2.67]

Viscosity_products = zeros(length(products))
Viscosity_products[4] = 31.5

Sulphur_3 =[	0.1]

BarrelToKT = zeros(length(crudes))
Sulphur_GO_stdev = zeros(length(crudes))
VaccuumResidue_nominal = zeros(length(crudes))
VaccuumResidue_stdev = zeros(length(crudes))

BarrelToKT = (GranularityOfBarrels/6.29)*(Crude_density/1000);
Sulphur_GO_stdev = 0.1*Sulphur_GO_nominal;
VaccuumResidue_nominal = Crude_yield[:,length(components)];
VaccuumResidue_stdev = 0.1*VaccuumResidue_nominal;


Crude_lower_bound =  (Barrel_lower_bound/GranularityOfBarrels)*BarrelToKT;
Crude_upper_bound =   (Barrel_upper_bound/GranularityOfBarrels)*BarrelToKT;


scenarios = 1:5
include("refinery1_data/5.jl")

Crude_yield_data = zeros(length(crudes), length(components), length(scenarios))
Desulphurisation_cost = zeros(length(crudes), length(scenarios))
Sulphur_2 = zeros(length(crudes), length(scenarios))

for h in scenarios
	for c in crudes
		for w in components
			if w < length(components)
				Crude_yield_data[c,w,h] = Crude_yield[c,w]/(1-Crude_yield[c,length(components)])*(1-VaccuumResidue_data[c,h])
			else
				Crude_yield_data[c,w,h] = VaccuumResidue_data[c,h]
			end
		end
		Desulphurisation_cost[c,h] = ((Sulphur_GO_data[c,h]*109.0909 + 365.4546)/1000)*(0.85/0.159)/GO_density[c]
		Sulphur_2[c,h] = Sulphur_GO_data[c,h]*0.005
	end
end






















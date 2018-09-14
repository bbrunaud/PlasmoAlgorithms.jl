include("input.jl")
function generate_model()
	
	m = Model(solver=BaronSolver(maxtime=1e4, epsr= 1e-5, CplexLibName = "/opt/ibm/ILOG/CPLEX_Studio127/cplex/bin/x86-64_linux/libcplex1270.so"))
	@variable(m, pickCrude[c in crudes], Bin)

	@variable(m, crudeQuantity[c in crudes]>=0)
	@variable(m, slack1[c in crudes, h in scenarios]>=0)
	@variable(m, slack2[c in crudes, h in scenarios]>=0)
	@variable(m, Reformer95_lower<=flow_Reformer95[h in scenarios]<=Reformer_capacity)
	@variable(m, 0<=flow_Reformer100[h in scenarios]<=Reformer_capacity - Reformer95_lower)
	@variable(m, 0<=flow_Cracker_Mogas[h in scenarios]<=Cracker_capacity)
	@variable(m, 0<=flow_Cracker_AGO[h in scenarios]<=Cracker_capacity)
	@variable(m, 0<=flow_Isomerisation[h in scenarios]<=CDU_capacity)
	@variable(m, 0<=flow_Desulphurisation_CGO[h in scenarios]<=Cracker_capacity)
	@variable(m, 0<=flow_LG_producing[h in scenarios]<=CDU_capacity)
	@variable(m, 0<=flow_LN_producing[h in scenarios]<=CDU_capacity)
	@variable(m, 0<=flow_HF_2[h in scenarios]<=Cracker_capacity)
	@variable(m, 0<=volume_PG98[h in scenarios]<=CDU_capacity/Density_PG98_input[1])
	@variable(m, 0<=volume_ES95[h in scenarios]<=CDU_capacity/Density_PG98_input[1])
	@variable(m, 0<=volume_HF[h in scenarios]<=CDU_capacity/GO_density[7])

	@variable(m, 0<=blin_CDU_LG[k in LG_out,h in scenarios]<=CDU_capacity)
	@variable(m, 0<=blin_Reformer95_LG[k in LG_out,h in scenarios]<=CDU_capacity)
	@variable(m, 0<=blin_Reformer100_LG[k in LG_out,h in scenarios]<=CDU_capacity)
	@variable(m, 0<=blin_Mogas_LG[k in LG_out,h in scenarios]<=CDU_capacity)
	@variable(m, 0<=blin_AGO_LG[k in LG_out,h in scenarios]<=CDU_capacity)
	@variable(m, 0<=blin_Cracker_Mogas[k in Cr_CGO,h in scenarios]<=Cracker_capacity)
	@variable(m, 0<=blin_Cracker_AGO[k in Cr_CGO,h in scenarios]<=Cracker_capacity)

	@variable(m, 0<=flow_Desulphurisation_1[c in crudes, h in scenarios]<=Desulphurisation_capacity)
	@variable(m, 0<=flow_AGO_1[c in crudes, h in scenarios]<=Crude_upper_bound[c])
	@variable(m, 0<=flow_AGO_2[c in crudes, h in scenarios]<=Desulphurisation_capacity)
	@variable(m, 0<=flow_HF_1[c in crudes, h in scenarios]<=Crude_upper_bound[c])
	@variable(m, 0<=flow_HF_3[c in crudes, h in scenarios]<=Crude_upper_bound[c])
	@variable(m, 0<=flow_Burn[k in Burn,h in scenarios]<=CDU_capacity)
	@variable(m, 0<=flow_PG98[k in PG98_in,h in scenarios]<=CDU_capacity)
	@variable(m, 0<=flow_ES95[k in PG98_in,h in scenarios]<=CDU_capacity)
	@variable(m, 0<=flow_AGO_3[k in AGO_in,h in scenarios]<=Cracker_capacity)
	setupperbound.(flow_AGO_3[1, :], CDU_capacity)
	@variable(m, 0<=flow_JPF[k in JPF_out,h in scenarios]<=CDU_capacity)
	@variable(m, 0<=flow_Import[p in products,h in scenarios]<=Import_upper[p])
	@variable(m, 0<=fraction_LG[k in LG_in,h in scenarios]<=1)
	@variable(m, 0<=fraction_CGO[k in Cr_mode,h in scenarios]<=1)








	@constraint(m, CDU_capacity_bound, 	sum(crudeQuantity[c]*BarrelToKT[c]/GranularityOfBarrels for c in crudes) <= CDU_capacity)


	@constraint(m, Crude_selection[c in crudes], crudeQuantity[c] >= pickCrude[c]*Barrel_lower_bound)


	@constraint(m, Crude_bound[c in crudes], 	crudeQuantity[c] <= pickCrude[c]*Barrel_upper_bound)


	@constraint(m, Desulphurisation_capacity_bound[h in scenarios], flow_Desulphurisation_CGO[h] + sum(flow_Desulphurisation_1[c,h] for c in crudes) <= Desulphurisation_capacity)


	@constraint(m, Mass_balance1[h in scenarios], Reformer_fraction[1,1]*flow_Reformer95[h] +
						Reformer_fraction[2,1]*flow_Reformer100[h] +
						Cracker_fraction[1,1]*flow_Cracker_Mogas[h] +
						Cracker_fraction[2,1]*flow_Cracker_AGO[h] +
						Isomerisation_fraction[1]*flow_Isomerisation[h] +
						Desulphurisation_fraction2[2]*flow_Desulphurisation_CGO[h] -
						flow_Burn[1,h] +
						sum(
							Crude_yield_data[c,1,h]*(crudeQuantity[c] + slack1[c,h] - slack2[c,h])*BarrelToKT[c]/GranularityOfBarrels +
							Desulphurisation_fraction[c,2]*flow_Desulphurisation_1[c,h]  for c in crudes
						) == 0)


	@constraint(m, Mass_balance2[h in scenarios], Reformer_fraction[1,2]*flow_Reformer95[h] +
						Reformer_fraction[2,2]*flow_Reformer100[h] +
						Cracker_fraction[1,2]*flow_Cracker_Mogas[h] +
						Cracker_fraction[2,2]*flow_Cracker_AGO[h] -
						flow_LG_producing[h] - flow_PG98[1,h] -
						flow_ES95[1,h] - flow_Burn[2,h] +
						sum(
							Crude_yield_data[c,2,h]*(crudeQuantity[c] + slack1[c,h] - slack2[c,h])*BarrelToKT[c]/GranularityOfBarrels for c in crudes
						) == 0)


	@constraint(m, Mass_balance3[h in scenarios], -flow_LN_producing[h] - flow_Burn[3,h] -
						flow_PG98[3,h] - flow_ES95[3,h] -
						flow_Isomerisation[h] - flow_JPF[1,h]*JPF_fraction[1,1] -
						flow_JPF[2,h]*JPF_fraction[1,2] +
						sum(
							Crude_yield_data[c,3,h]*(crudeQuantity[c] + slack1[c,h] - slack2[c,h])*BarrelToKT[c]/GranularityOfBarrels for c in crudes
						) == 0)


	@constraint(m, Mass_balance4[h in scenarios], -flow_JPF[1,h]*JPF_fraction[2,1] -
						flow_JPF[2,h]*JPF_fraction[2,2] -
						flow_Reformer95[h] - flow_Reformer100[h] +
						sum(
							Crude_yield_data[c,4,h]*(crudeQuantity[c] + slack1[c,h] - slack2[c,h])*BarrelToKT[c]/GranularityOfBarrels for c in crudes
						) == 0)


	@constraint(m, Mass_balance5[h in scenarios], -flow_JPF[1,h]*JPF_fraction[3,1] -
						flow_JPF[2,h]*JPF_fraction[3,2] -
						flow_AGO_3[1,h] +
						sum(
							Crude_yield_data[c,5,h]*(crudeQuantity[c] + slack1[c,h] - slack2[c,h])*BarrelToKT[c]/GranularityOfBarrels for c in crudes
						) == 0)


	@constraint(m, Mass_balance7[h in scenarios], -flow_Cracker_Mogas[h] - flow_Cracker_AGO[h] +
						sum(
							Crude_yield_data[c,7,h]*(crudeQuantity[c] + slack1[c,h] - slack2[c,h])*BarrelToKT[c]/GranularityOfBarrels for c in crudes
						) == 0)


	@constraint(m, GO_balance[c in crudes,h in scenarios], -flow_AGO_1[c,h] - flow_Desulphurisation_1[c,h] - flow_HF_3[c,h] +
						Crude_yield_data[c,6,h]*(crudeQuantity[c] + slack1[c,h] - slack2[c,h])*BarrelToKT[c]/GranularityOfBarrels == 0)


	@constraint(m, VR_balance[c in crudes,h in scenarios], 	Crude_yield_data[c,8,h]*(crudeQuantity[c] + slack1[c,h] - slack2[c,h])*BarrelToKT[c]/GranularityOfBarrels == flow_HF_1[c,h])


	@constraint(m, Desulphurisation_balance[c in crudes,h in scenarios], Desulphurisation_fraction[c,1]*flow_Desulphurisation_1[c,h] == flow_AGO_2[c,h])


	@constraint(m, Reformer95_balance[h in scenarios], 	flow_Reformer95[h]*Reformer_fraction[1,3] +
								flow_Reformer100[h]*Reformer_fraction[2,3] ==
								flow_PG98[4,h] + flow_ES95[4,h])


	@constraint(m, Reformer100_balance[h in scenarios], 	flow_Reformer95[h]*Reformer_fraction[1,4] +
								flow_Reformer100[h]*Reformer_fraction[2,4] ==
								flow_PG98[5,h] + flow_ES95[5,h])


	@constraint(m, Isomerisation_balance[h in scenarios], flow_Isomerisation[h]*Isomerisation_fraction[2] ==
								flow_PG98[2,h] + flow_ES95[2,h])


	@constraint(m, CN_balance[h in scenarios], 	flow_Cracker_Mogas[h]*Cracker_fraction[1,3] +
						flow_Cracker_AGO[h]*Cracker_fraction[2,3] ==
						flow_PG98[6,h] + flow_ES95[6,h])


	@constraint(m, CGO_balance[h in scenarios], 	flow_Cracker_Mogas[h]*Cracker_fraction[1,4] +
						flow_Cracker_AGO[h]*Cracker_fraction[2,4] ==
						flow_Desulphurisation_CGO[h] + flow_HF_2[h] + flow_AGO_3[2,h])


	@constraint(m, Desulphurisation_CGO_balance[h in scenarios], Desulphurisation_fraction2[1]*flow_Desulphurisation_CGO[h] == flow_AGO_3[3,h])


	@constraint(m, Demand_constraint1[h in scenarios], flow_Import[1,h] + sum(flow_PG98[k,h] for k in PG98_in) >= Demand_quantity[1])


	@constraint(m, Demand_constraint2[h in scenarios], flow_Import[2,h] + sum(flow_ES95[k,h] for k in PG98_in) >= Demand_quantity[2])


	@constraint(m, Demand_constraint3[h in scenarios], flow_Import[3,h] + sum(flow_JPF[k,h] for k in JPF_out) >= Demand_quantity[3])


	@constraint(m, Demand_constraint4[h in scenarios], flow_Import[4,h] + sum(flow_AGO_3[k,h] for k in AGO_in) +
								sum(flow_AGO_1[c,h] + flow_AGO_2[c,h] for c in crudes) >= Demand_quantity[4])


	@constraint(m, Demand_constraint5[h in scenarios], flow_Import[5,h] + flow_HF_2[h] + 
								sum(flow_HF_1[c,h] + flow_HF_3[c,h] for c in crudes) >= Demand_quantity[5])


	@constraint(m, Demand_constraint6[h in scenarios], flow_Import[6,h] + flow_LG_producing[h] >= Demand_quantity[6])


	@constraint(m, Demand_constraint7[h in scenarios], flow_Import[7,h] + flow_LN_producing[h] >= Demand_quantity[7])


	@constraint(m, PG98_volume_def[h in scenarios],	flow_Import[1,h]/Density_products[1] +
							sum(flow_PG98[k,h]/Density_PG98_input[k] for k in PG98_in) == volume_PG98[h])


	@constraint(m, ES95_volume_def[h in scenarios],	flow_Import[2,h]/Density_products[2] +
							sum(flow_ES95[k,h]/Density_PG98_input[k] for k in PG98_in) == volume_ES95[h])


	@constraint(m, Butane95_constraint[h in scenarios],	flow_ES95[1,h]/Density_PG98_input[1] +
								0.03*flow_Import[2,h]/Density_products[2] <= 0.05*volume_ES95[h])


	@constraint(m, Butane98_constraint[h in scenarios],	flow_PG98[1,h]/Density_PG98_input[1] +
								0.03*flow_Import[1,h]/Density_products[2] <= 0.05*volume_PG98[h])


	@NLconstraint(m, blincon_CDU_LG1[h in scenarios], blin_CDU_LG[1,h] == fraction_LG[1,h]*flow_ES95[1,h])


	@NLconstraint(m, blincon_CDU_LG2[h in scenarios], blin_CDU_LG[2,h] == fraction_LG[1,h]*flow_PG98[1,h])


	@NLconstraint(m, blincon_CDU_LG3[h in scenarios], blin_CDU_LG[3,h] == fraction_LG[1,h]*flow_Burn[2,h])


	@NLconstraint(m, blincon_CDU_LG4[h in scenarios], blin_CDU_LG[4,h] == fraction_LG[1,h]*flow_LG_producing[h])


	@NLconstraint(m, blincon_Reformer95_LG1[h in scenarios], blin_Reformer95_LG[1,h] == fraction_LG[2,h]*flow_ES95[1,h])


	@NLconstraint(m, blincon_Reformer95_LG2[h in scenarios], blin_Reformer95_LG[2,h] == fraction_LG[2,h]*flow_PG98[1,h])


	@NLconstraint(m, blincon_Reformer95_LG3[h in scenarios], blin_Reformer95_LG[3,h] == fraction_LG[2,h]*flow_Burn[2,h])


	@NLconstraint(m, blincon_Reformer95_LG4[h in scenarios], blin_Reformer95_LG[4,h] == fraction_LG[2,h]*flow_LG_producing[h])


	@NLconstraint(m, blincon_Reformer100_LG1[h in scenarios], blin_Reformer100_LG[1,h] == fraction_LG[3,h]*flow_ES95[1,h])


	@NLconstraint(m, blincon_Reformer100_LG2[h in scenarios], blin_Reformer100_LG[2,h] == fraction_LG[3,h]*flow_PG98[1,h])


	@NLconstraint(m, blincon_Reformer100_LG3[h in scenarios], blin_Reformer100_LG[3,h] == fraction_LG[3,h]*flow_Burn[2,h])


	@NLconstraint(m, blincon_Reformer100_LG4[h in scenarios], blin_Reformer100_LG[4,h] == fraction_LG[3,h]*flow_LG_producing[h])


	@NLconstraint(m, blincon_Mogas_LG1[h in scenarios], blin_Mogas_LG[1,h] == fraction_LG[4,h]*flow_ES95[1,h])


	@NLconstraint(m, blincon_Mogas_LG2[h in scenarios], blin_Mogas_LG[2,h] == fraction_LG[4,h]*flow_PG98[1,h])


	@NLconstraint(m, blincon_Mogas_LG3[h in scenarios], blin_Mogas_LG[3,h] == fraction_LG[4,h]*flow_Burn[2,h])


	@NLconstraint(m, blincon_Mogas_LG4[h in scenarios], blin_Mogas_LG[4,h] == fraction_LG[4,h]*flow_LG_producing[h])


	@NLconstraint(m, blincon_AGO_LG1[h in scenarios], blin_AGO_LG[1,h] == fraction_LG[5,h]*flow_ES95[1,h])


	@NLconstraint(m, blincon_AGO_LG2[h in scenarios], blin_AGO_LG[2,h] == fraction_LG[5,h]*flow_PG98[1,h])


	@NLconstraint(m, blincon_AGO_LG3[h in scenarios], blin_AGO_LG[3,h] == fraction_LG[5,h]*flow_Burn[2,h])


	@NLconstraint(m, blincon_AGO_LG4[h in scenarios], blin_AGO_LG[4,h] == fraction_LG[5,h]*flow_LG_producing[h])


	@constraint(m, LG_balance[h in scenarios], 	sum(blin_CDU_LG[k,h] for k in LG_out) ==
						sum(
							Crude_yield_data[c,2,h]*(crudeQuantity[c] + slack1[c,h] - slack2[c,h])*BarrelToKT[c]/GranularityOfBarrels
						for c in crudes))


	@constraint(m, Reformer95_LG_balance[h in scenarios], flow_Reformer95[h]*Reformer_fraction[1,2] ==
								sum(blin_Reformer95_LG[k,h] for k in LG_out))


	@constraint(m, Reformer100_LG_balance[h in scenarios], 	flow_Reformer100[h]*Reformer_fraction[2,2] ==
									sum(blin_Reformer100_LG[k,h] for k in LG_out))


	@constraint(m, Cracker_Mogas_LG_balance[h in scenarios], 	flow_Cracker_Mogas[h]*Cracker_fraction[1,2] ==
									sum(blin_Mogas_LG[k,h] for k in LG_out))


	@constraint(m, Cracker_AGO_LG_balance[h in scenarios], 	flow_Cracker_AGO[h]*Cracker_fraction[2,2] ==
									sum(blin_AGO_LG[k,h] for k in LG_out))


	@constraint(m, pq_ES95_constraint[h in scenarios], 	blin_CDU_LG[1,h] + blin_Reformer95_LG[1,h] +
								blin_Reformer100_LG[1,h] + blin_Mogas_LG[1,h] +
								blin_AGO_LG[1,h] == flow_ES95[1,h])


	@constraint(m, pq_PG98_constraint[h in scenarios], 	blin_CDU_LG[2,h] + blin_Reformer95_LG[2,h] +
								blin_Reformer100_LG[2,h] + blin_Mogas_LG[2,h] +
								blin_AGO_LG[2,h] == flow_PG98[1,h])


	@constraint(m, pq_burn_constraint[h in scenarios], 	blin_CDU_LG[3,h] + blin_Reformer95_LG[3,h] +
								blin_Reformer100_LG[3,h] + blin_Mogas_LG[3,h] +
								blin_AGO_LG[3,h] == flow_Burn[2,h])


	@constraint(m, pq_demand_constraint[h in scenarios], 	blin_CDU_LG[4,h] + blin_Reformer95_LG[4,h] +
								blin_Reformer100_LG[4,h] + blin_Mogas_LG[4,h] +
								blin_AGO_LG[4,h] == flow_LG_producing[h])


	@constraint(m, LG_split_balance[h in scenarios], sum(fraction_LG[k,h] for k in LG_in) == 1)


	@constraint(m, VP_ES95_lower[h in scenarios], -0.45*volume_ES95[h] + flow_Import[2,h]*Product_VP[2]/Density_products[2] +
						sum(VP[k]*flow_ES95[k,h]/Density_PG98_input[k] for k in PG98_in) +
						LG_parameters[1,1]*blin_CDU_LG[1,h]/Density_PG98_input[1] +
						LG_parameters[1,2]*blin_Reformer95_LG[1,h]/Density_PG98_input[1] +
						LG_parameters[1,3]*blin_Reformer100_LG[1,h]/Density_PG98_input[1] +
						LG_parameters[1,4]*blin_Mogas_LG[1,h]/Density_PG98_input[1] +
						LG_parameters[1,5]*blin_AGO_LG[1,h]/Density_PG98_input[1] >= 0)


	@constraint(m, VP_ES95_upper[h in scenarios], -0.80*volume_ES95[h] + flow_Import[2,h]*Product_VP[2]/Density_products[2] +
						sum(VP[k]*flow_ES95[k,h]/Density_PG98_input[k] for k in PG98_in) +
						LG_parameters[1,1]*blin_CDU_LG[1,h]/Density_PG98_input[1] +
						LG_parameters[1,2]*blin_Reformer95_LG[1,h]/Density_PG98_input[1] +
						LG_parameters[1,3]*blin_Reformer100_LG[1,h]/Density_PG98_input[1] +
						LG_parameters[1,4]*blin_Mogas_LG[1,h]/Density_PG98_input[1] +
						LG_parameters[1,5]*blin_AGO_LG[1,h]/Density_PG98_input[1] <= 0)


	@constraint(m, VP_PG98_lower[h in scenarios], -0.50*volume_PG98[h] + flow_Import[1,h]*Product_VP[1]/Density_products[1] +
						sum(VP[k]*flow_PG98[k,h]/Density_PG98_input[k] for k in PG98_in) +
						LG_parameters[1,1]*blin_CDU_LG[2,h]/Density_PG98_input[1] +
						LG_parameters[1,2]*blin_Reformer95_LG[2,h]/Density_PG98_input[1] +
						LG_parameters[1,3]*blin_Reformer100_LG[2,h]/Density_PG98_input[1] +
						LG_parameters[1,4]*blin_Mogas_LG[2,h]/Density_PG98_input[1] +
						LG_parameters[1,5]*blin_AGO_LG[2,h]/Density_PG98_input[1] >= 0)


	@constraint(m, VP_PG98_upper[h in scenarios], -0.86*volume_PG98[h] + flow_Import[1,h]*Product_VP[1]/Density_products[1] +
						sum(VP[k]*flow_PG98[k,h]/Density_PG98_input[k] for k in PG98_in) +
						LG_parameters[1,1]*blin_CDU_LG[2,h]/Density_PG98_input[1] +
						LG_parameters[1,2]*blin_Reformer95_LG[2,h]/Density_PG98_input[1] +
						LG_parameters[1,3]*blin_Reformer100_LG[2,h]/Density_PG98_input[1] +
						LG_parameters[1,4]*blin_Mogas_LG[2,h]/Density_PG98_input[1] +
						LG_parameters[1,5]*blin_AGO_LG[2,h]/Density_PG98_input[1] <= 0)


	@constraint(m, RON_PG98[h in scenarios], 	-98*volume_PG98[h] + flow_Import[1,h]*Product_RON[1]/Density_products[1] +
					sum( RON[k]*flow_PG98[k,h]/Density_PG98_input[k] for k in PG98_in) +
					LG_parameters[2,1]*blin_CDU_LG[2,h]/Density_PG98_input[1] +
					LG_parameters[2,2]*blin_Reformer95_LG[2,h]/Density_PG98_input[1] +
					LG_parameters[2,3]*blin_Reformer100_LG[2,h]/Density_PG98_input[1] +
					LG_parameters[2,4]*blin_Mogas_LG[2,h]/Density_PG98_input[1] +
					LG_parameters[2,5]*blin_AGO_LG[2,h]/Density_PG98_input[1] >= 0)


	@constraint(m, RON_ES95[h in scenarios], 	-95*volume_ES95[h] + flow_Import[2,h]*Product_RON[2]/Density_products[2] +
					sum(RON[k]*flow_ES95[k,h]/Density_PG98_input[k] for k in PG98_in) +
					LG_parameters[2,1]*blin_CDU_LG[1,h]/Density_PG98_input[1] +
					LG_parameters[2,2]*blin_Reformer95_LG[1,h]/Density_PG98_input[1] +
					LG_parameters[2,3]*blin_Reformer100_LG[1,h]/Density_PG98_input[1] +
					LG_parameters[2,4]*blin_Mogas_LG[1,h]/Density_PG98_input[1] +
					LG_parameters[2,5]*blin_AGO_LG[1,h]/Density_PG98_input[1] >= 0)


	@constraint(m, Sensitivity_PG98[h in scenarios], 	-10*volume_PG98[h] + flow_Import[1,h]*(Product_RON[1] - Product_MON[1])/Density_products[1] +
							sum((RON[k] - MON[k])*flow_PG98[k,h]/Density_PG98_input[k] for k in PG98_in) +
							(LG_parameters[2,1] - LG_parameters[3,1])*blin_CDU_LG[2,h]/Density_PG98_input[1] +
							(LG_parameters[2,2] - LG_parameters[3,2])*blin_Reformer95_LG[2,h]/Density_PG98_input[1] +
							(LG_parameters[2,3] - LG_parameters[3,3])*blin_Reformer100_LG[2,h]/Density_PG98_input[1] +
							(LG_parameters[2,4] - LG_parameters[3,4])*blin_Mogas_LG[2,h]/Density_PG98_input[1] +
							(LG_parameters[2,5] - LG_parameters[3,5])*blin_AGO_LG[2,h]/Density_PG98_input[1] <= 0)


	@constraint(m, Sensitivity_ES95[h in scenarios], 	-10*volume_ES95[h] + flow_Import[2,h]*(Product_RON[2] - Product_MON[2])/Density_products[2] +
							sum((RON[k] - MON[k])*flow_ES95[k,h]/Density_PG98_input[k] for k in PG98_in) +
							(LG_parameters[2,1] - LG_parameters[3,1])*blin_CDU_LG[1,h]/Density_PG98_input[1] +
							(LG_parameters[2,2] - LG_parameters[3,2])*blin_Reformer95_LG[1,h]/Density_PG98_input[1] +
							(LG_parameters[2,3] - LG_parameters[3,3])*blin_Reformer100_LG[1,h]/Density_PG98_input[1] +
							(LG_parameters[2,4] - LG_parameters[3,4])*blin_Mogas_LG[1,h]/Density_PG98_input[1] +
							(LG_parameters[2,5] - LG_parameters[3,5])*blin_AGO_LG[1,h]/Density_PG98_input[1] <= 0)


	@NLconstraint(m, blincon_Cracker_Mogas1[h in scenarios], blin_Cracker_Mogas[1,h] == fraction_CGO[1,h]*flow_AGO_3[2,h])


	@NLconstraint(m, blincon_Cracker_Mogas2[h in scenarios], blin_Cracker_Mogas[2,h] == fraction_CGO[1,h]*flow_HF_2[h])


	@NLconstraint(m, blincon_Cracker_Mogas3[h in scenarios], blin_Cracker_Mogas[3,h] == fraction_CGO[1,h]*flow_Desulphurisation_CGO[h])


	@NLconstraint(m, blincon_Cracker_AGO1[h in scenarios], blin_Cracker_AGO[1,h] == fraction_CGO[2,h]*flow_AGO_3[2,h])


	@NLconstraint(m, blincon_Cracker_AGO2[h in scenarios], blin_Cracker_AGO[2,h] == fraction_CGO[2,h]*flow_HF_2[h])


	@NLconstraint(m, blincon_Cracker_AGO3[h in scenarios], blin_Cracker_AGO[3,h] == fraction_CGO[2,h]*flow_Desulphurisation_CGO[h])


	@constraint(m, Cracker_Mogas_CGO_balance[h in scenarios], blin_Cracker_Mogas[1,h] + blin_Cracker_Mogas[2,h] +
									blin_Cracker_Mogas[3,h] == flow_Cracker_Mogas[h]*Cracker_fraction[1,4])


	@constraint(m, Cracker_AGO_CGO_balance[h in scenarios], 	blin_Cracker_AGO[1,h] + blin_Cracker_AGO[2,h] +
									blin_Cracker_AGO[3,h] == flow_Cracker_AGO[h]*Cracker_fraction[2,4])


	@constraint(m, CGO_split_balance[h in scenarios], sum(fraction_CGO[k,h] for k in Cr_mode) == 1)


	@constraint(m, pq_AGO_constraint[h in scenarios], blin_Cracker_Mogas[1,h] + blin_Cracker_AGO[1,h] == flow_AGO_3[2,h])


	@constraint(m, pq_HF_constraint[h in scenarios], blin_Cracker_Mogas[2,h] + blin_Cracker_AGO[2,h] == flow_HF_2[h])


	@constraint(m, pq_Desulphurisation_constraint[h in scenarios], blin_Cracker_Mogas[3,h] + blin_Cracker_AGO[3,h] == flow_Desulphurisation_CGO[h])


	@constraint(m, HF_volume_def[h in scenarios], -volume_HF[h] + flow_Import[5,h]/Density_products[5] +
						flow_HF_2[h]/CGO_density +
						sum(flow_HF_1[c,h]/HFO_density[c] + flow_HF_3[c,h]/GO_density[c] for c in crudes) == 0)


	@constraint(m, HF_viscosity_lower[h in scenarios], 	flow_Import[5,h]*Viscosity_products[5]/Density_products[5] +
								sum(
									flow_HF_1[c,h]*Viscosity_HF1[c]/HFO_density[c] +
									flow_HF_3[c,h]*Viscosity_HF3[c]/GO_density[c] for c in crudes
								) +
								(blin_Cracker_Mogas[2,h]*Mogas_viscosity + blin_Cracker_AGO[2,h]*AGO_viscosity)/CGO_density -
								30*volume_HF[h] >= 0)


	@constraint(m, HF_viscosity_upper[h in scenarios], 	flow_Import[5,h]*Viscosity_products[5]/Density_products[5] +
								sum(
									flow_HF_1[c,h]*Viscosity_HF1[c]/HFO_density[c] +
									flow_HF_3[c,h]*Viscosity_HF3[c]/GO_density[c] for c in crudes
								) +
								(blin_Cracker_Mogas[2,h]*Mogas_viscosity + blin_Cracker_AGO[2,h]*AGO_viscosity)/CGO_density -
								33*volume_HF[h] <= 0)


	@constraint(m, AGO_sulphur_balance[h in scenarios], 	flow_Import[4,h]*Product_Sulphur[4] - Sulphur_spec*flow_Import[4,h] +
								sum(
									(Sulphur_GO_data[c,h] - Sulphur_spec)*flow_AGO_1[c,h] +
									(Sulphur_2[c,h] - Sulphur_spec)*flow_AGO_2[c,h] for c in crudes
								) +
								flow_AGO_3[1,h]*(Sulphur_3[1] - Sulphur_spec) +
								blin_Cracker_AGO[1,h]*(AGO_Sulphur - Sulphur_spec) +
								blin_Cracker_Mogas[1,h]*(Mogas_Sulphur - Sulphur_spec) +
								blin_Cracker_AGO[3,h]*AGO_Sulphur*0.005 +
								blin_Cracker_Mogas[3,h]*Mogas_Sulphur*0.005 -
								Sulphur_spec*flow_AGO_3[3,h] <= 0)


	@constraint(m, Refinery_Fuel[h in scenarios], 1.3*flow_Burn[1,h] + 1.2*flow_Burn[2,h] + 1.1*flow_Burn[3,h] -
						flow_Reformer95[h]*Reformer_fraction[1,5] -
						flow_Reformer100[h]*Reformer_fraction[2,5] -
						flow_Cracker_Mogas[h]*Cracker_fraction[1,5] -
						flow_Cracker_AGO[h]*Cracker_fraction[2,5] -
						flow_Isomerisation[h]*Isomerisation_fraction[3] -
						flow_Desulphurisation_CGO[h]*Desulphurisation_fraction2[3] - 15.2 -
						sum(
							0.018*(crudeQuantity[c] + slack1[c,h] - slack2[c,h])*BarrelToKT[c]/GranularityOfBarrels +
							flow_Desulphurisation_1[c,h]*Desulphurisation_fraction[c,3] for c in crudes
						) >= 0)


	@constraint(m, Cracker_capacity_bound[h in scenarios], flow_Cracker_Mogas[h] + flow_Cracker_AGO[h] <= Cracker_capacity)


	@constraint(m, Reformer_capacity_bound[h in scenarios], flow_Reformer95[h] + flow_Reformer100[h] <= Reformer_capacity)



   	@objective(m, Min, sum(prob[h]*(
								Cracker_Mogas_cost*flow_Cracker_Mogas[h] +
								Cracker_AGO_cost*flow_Cracker_AGO[h] +
								Reformer95_cost*flow_Reformer95[h] +
								Reformer100_cost*flow_Reformer100[h] +
								Isomerisation_cost*flow_Isomerisation[h] +
								Desulphurisation_CGO_cost*flow_Desulphurisation_CGO[h] -
								LG_sale*flow_LG_producing[h] -
								LN_sale*flow_LN_producing[h] -
								HF_sale*flow_HF_2[h] +
								sum(
									Desulphurisation_cost[c,h]*flow_Desulphurisation_1[c,h] -
									AGO_sale*flow_AGO_1[c,h] -
									AGO_sale*flow_AGO_2[c,h] -
									HF_sale*flow_HF_1[c,h] -
									HF_sale*flow_HF_3[c,h] +
									((crudeQuantity[c] + 100 * slack1[c,h] + 100*slack2[c,h])/1000)*(Crude_price[c]+1) for c in crudes
								) -
								sum(
									PG98_sale*flow_PG98[k,h] +
									ES95_sale*flow_ES95[k,h] for k in PG98_in
								) -
								sum(
									JET_sale*flow_JPF[k,h] for k in JPF_out
								) -
								sum(
									AGO_sale*flow_AGO_3[k,h] for k in AGO_in
								)
							) for h in scenarios))

   	return m
end

model = generate_model()
solve(model)
function generate_master()
	
	m = Model(solver=CplexSolver(CPX_PARAM_SCRIND=0))
	@variable(m, pickCrude[c in crudes], Bin)

	@variable(m, crudeQuantity[c in crudes]>=0)
	







	@constraint(m, CDU_capacity_bound, 	sum(crudeQuantity[c]*BarrelToKT[c]/GranularityOfBarrels for c in crudes) <= CDU_capacity)


	@constraint(m, Crude_selection[c in crudes], crudeQuantity[c] >= pickCrude[c]*Barrel_lower_bound)


	@constraint(m, Crude_bound[c in crudes], 	crudeQuantity[c] <= pickCrude[c]*Barrel_upper_bound)


	

   	@objective(m, Min, sum((crudeQuantity[c]/1000)*(Crude_price[c]+1) for c in crudes))

   	return m
end


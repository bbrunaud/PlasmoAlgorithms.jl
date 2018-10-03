function generate_master()
	
	m = Model(solver=CplexSolver(CPX_PARAM_SCRIND=0))
	@variable(m, pickCrude[c in crudes], Bin)

	@variable(m, 0<=crudeQuantity[c in crudes]<=floor(Crude_upper_bound[c]*1000)/1000)
	







	@constraint(m, CDU_capacity_bound, 	sum(crudeQuantity[c] for c in crudes) <= CDU_capacity)


	@constraint(m, Crude_selection[c in crudes], crudeQuantity[c] >= pickCrude[c]*ceil(Crude_lower_bound[c]*1000)/1000)


	@constraint(m, Crude_bound[c in crudes], 	crudeQuantity[c] <= pickCrude[c]*floor(Crude_upper_bound[c]*1000)/1000)


	

   	@objective(m, Min, sum((crudeQuantity[c]/1000)/BarrelToKT[c]*GranularityOfBarrels*(Crude_price[c]+1) for c in crudes))

   	return m
end


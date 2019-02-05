function generate_sub(prob = prob, Return = Return, scenario = 1, is_laststage= is_laststage)
	Return_s = Return[:, scenario]
	sub = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0))
	@variable(sub, AmountInvestedPreviousStage[i in Investments]>=0)
	@variable(sub, AmountInvestedCurrentStage[i in Investments]>=0)

	@variable(sub, InvestmentType[i in Investments], Bin)

	@constraint(sub, sum(Return_s[i] * AmountInvestedPreviousStage[i] for i in Investments) == sum(AmountInvestedCurrentStage[i] for i in Investments))
	@constraint(sub, sum(InvestmentType[i] for i in Investments) == 1)
	@constraint(sub, c[i in Investments], AmountInvestedCurrentStage[i] <= 120000 * InvestmentType[i])

	if is_laststage
		@variable(sub, TargetDeficit[s in base_scenarios]>=0)
		@variable(sub, TargetSurplus[s in base_scenarios]>=0)
		@variable(sub, FinalWealth[s in base_scenarios])
		@constraint(sub, c1[s in base_scenarios], sum(AmountInvestedCurrentStage[i] * Return[i, s] for i in Investments) - TargetSurplus[s] + TargetDeficit[s] - TargetCapital == 0.0 )		
		@constraint(sub, c2[s in base_scenarios], FinalWealth[s] - ExcessPerUnitBenefit * TargetSurplus[s]  + ShortagePerUnitCost *  TargetDeficit[s] == 0.0 )
		@objective(sub, Min, -prob * sum(FinalWealth[s]/ length(base_scenarios) for s in base_scenarios))
	else
		@objective(sub, Min, 0)

	end
	return sub 
end
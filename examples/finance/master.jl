function generate_master()
	mp = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0))
	@variable(mp, AmountInvestedPreviousStage[i in Investments]>=0)
	@variable(mp, AmountInvestedCurrentStage[i in Investments]>=0)
	@variable(mp, InvestmentType[i in Investments], Bin)
	@constraint(mp, sum(InvestmentType[i] for i in Investments) == 1)
	@constraint(mp, c[i in Investments], AmountInvestedCurrentStage[i] <= 120000 * InvestmentType[i])
	@constraint(mp, sum(AmountInvested[i] for i in Investments) == InitialWealth)
	return mp
end
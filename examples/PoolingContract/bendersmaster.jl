function generate_bendersmaster()
	m = Model(solver=CplexSolver())

	@variable(m, gamma_intlt[i in feeds], Bin)
	@variable(m, gamma_pool[l in pools], Bin)
	@variable(m, S[l in pools]>=0)
	@variable(m, A[i in feeds]>=0)

	@constraint(m, f1[i in feeds], AL[i]*gamma_intlt[i]<= A[i])
	@constraint(m, f2[i in feeds], A[i] <= AU[i]*gamma_intlt[i])
	@constraint(m, f3[l in pools], SL[l]* gamma_pool[l]<=S[l])
	@constraint(m, f4[l in pools], S[l]<= SU[l]* gamma_pool[l])

	@objective(m, Min,  sum(c_fixed_inlt[i] * gamma_intlt[i]+c_variable_inlt[i]*A[i] for i in feeds) + sum(c_fixed_pool[l] * gamma_pool[l] + c_variable_pool[l]*S[l] for l in pools) )

	return m
end
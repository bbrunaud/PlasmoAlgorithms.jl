function crosssolve(bgraph::Plasmo.PlasmoGraph, lgraph::Plasmo.PlasmoGraph; 
	heuristic=:lagfirst, 
	max_iterations_lag::Int64=60, 
	max_iterations_benders::Int64=100,
	benders_cuts=[:LP],
	ϵ=1e-5,
	benders_UBupdatefrequency=1,
	benders_timelimit=3600,
	benders_verbose=false, 
	benders_LB=NaN, 
	is_nonconvex=false,
	lagrangeheuristic=:nearest_scenario, # function to calculate the upper bound
	lag_initialmultipliers=:zero, # :relaxation for LP relaxation
	lag_δ = 0.5, # Factor to shrink step when subgradient stuck
	lag_maxnoimprove = 1,
	lag_cpbound=1e6)# Amount of iterations that no improvement is allowed before shrinking step

	if heuristic == :lagfirst
		lagrangesolve(lgraph, max_iterations=max_iterations_lag, lagrangeheuristic=lagrangeheuristic,  maxnoimprove = 1)
		bendersolve(bgraph, cuts=benders_cuts, max_iterations=1, is_nonconvex=is_nonconvex)
		num_iter = length(getnodes(lgraph)[1].attributes[:Zsl])
		add_lagrangean_cuts(bgraph, lgraph, cut_indices=1:num_iter)
		bendersolve(bgraph, cuts=benders_cuts, max_iterations=1, is_nonconvex=is_nonconvex)
	end

end

function add_lagrangean_cuts(bgraph::Plasmo.PlasmoGraph, lgraph::Plasmo.PlasmoGraph; cut_indices =1:10)
	rootnode = bgraph.attributes[:roots][1]
	model = getmodel(rootnode)
	for node in values(getnodes(lgraph))
		Zsl = node.attributes[:Zsl][cut_indices]
		μ = node.attributes[:μ][cut_indices]
		index = rootnode.attributes[:scenario_to_θ][node.attributes[:scenario]]
		for i in 1:length(Zsl)
			rhs = Zsl[i]
			coeff = μ[i]
			for varname in keys(coeff)
				rhs -= coeff[varname] * rootnode.attributes[:varname_to_var][varname]
			end
			θ=getindex(model, :θ)
			@constraint(model, θ[index] >= rhs-node.attributes[:prob] *rootnode.attributes[:preobj].aff)
		end
	end

end
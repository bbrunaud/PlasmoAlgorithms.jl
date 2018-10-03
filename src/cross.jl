function crosssolve(bgraph::Plasmo.PlasmoGraph, lgraph::Plasmo.PlasmoGraph; 
	heuristic=:lagfirst, 
	max_iterations_lag::Int64=60, 
	max_iterations_benders::Int64=100,
	benders_cuts=[:LP],
	ϵ=1e-5,
	rel_gap=1e-4,
	benders_UBupdatefrequency=1,
	benders_timelimit=3600,
	benders_verbose=false, 
	benders_LB=NaN, 
	is_nonconvex=false,
	lag_α=2,
	lagrangeheuristic=nearest_scenario, # function to calculate the upper bound
	lag_initialmultipliers=:zero, # :relaxation for LP relaxation
	lag_δ = 0.5, # Factor to shrink step when subgradient stuck
	lag_maxnoimprove = 1,
	lag_cpbound=1e6,
	user_LB=-Inf,
	user_UB=+Inf)# Amount of iterations that no improvement is allowed before shrinking step
	results = Dict()
	results[:benders_time] = 0
	results[:lagrangean_time] = 0
	if heuristic == :lagfirst
		solution = lagrangesolve(lgraph, max_iterations=max_iterations_lag, ϵ=ϵ, rel_gap=rel_gap, initialmultipliers=lag_initialmultipliers, lagrangeheuristic=lagrangeheuristic,  maxnoimprove = lag_maxnoimprove, δ = lag_δ, α = lag_α, cpbound = lag_cpbound, user_LB=user_LB,user_UB=user_UB )
		println(solution.termination)
		results[:lagrangean_time] += solution.clocktime[end]
		if !lgraph.attributes[:is_infeasible]
			if !haskey(bgraph.attributes,:preprocessed)
				solution = bendersolve(bgraph, cuts=benders_cuts, max_iterations=1, is_nonconvex=is_nonconvex, timelimit=benders_timelimit, LB=benders_LB, ϵ=ϵ, rel_gap=rel_gap, UBupdatefrequency=benders_UBupdatefrequency, verbose=benders_verbose)
				results[:benders_time] += solution.clocktime[end]
				num_iter = length(getnodes(lgraph)[1].attributes[:Zsl])
				println("---------adding Lagrangean cuts to Benders master problem---------")
				add_lagrangean_cuts(bgraph, lgraph, cut_indices=1:num_iter)
				if max_iterations_benders < 2
					max_iterations_benders =2
				end
				solution = bendersolve(bgraph, cuts=benders_cuts, max_iterations=max_iterations_benders-1, is_nonconvex=is_nonconvex, timelimit=benders_timelimit, LB=benders_LB, ϵ=ϵ, rel_gap=rel_gap, UBupdatefrequency=benders_UBupdatefrequency, verbose=benders_verbose)
				results[:benders_time] += solution.clocktime[end]
			else
				num_iter = length(getnodes(lgraph)[1].attributes[:Zsl])
				println("---------adding Lagrangean cuts to Benders master problem---------")
				add_lagrangean_cuts(bgraph, lgraph, cut_indices=1:num_iter)
				solution = bendersolve(bgraph, cuts=benders_cuts, max_iterations=max_iterations_benders, is_nonconvex=is_nonconvex, timelimit=benders_timelimit, LB=benders_LB, ϵ=ϵ, rel_gap=rel_gap, UBupdatefrequency=benders_UBupdatefrequency, verbose=benders_verbose)
				results[:benders_time] += solution.clocktime[end]
			end
		else
			results[:LB] = +Inf 
			results[:UB] = 1e8
			results[:best_avg_x] = []
			results[:best_feasible_x] = []
			return results
		end
		results[:LB] = bgraph.attributes[:LB]
		results[:best_avg_x] = lgraph.attributes[:best_avg_x]
		if lgraph.attributes[:UB] < bgraph.attributes[:UB]
			results[:UB] = lgraph.attributes[:UB]
			results[:best_feasible_x] = lgraph.attributes[:best_feasible_x]
			results[:best_solution_source] = :lgraph
		else 
			results[:UB] = bgraph.attributes[:UB]
			results[:best_feasible_x] = bgraph.attributes[:best_feasible_x]
			results[:best_solution_source] = :bgraph
		end		
	end

	if heuristic == :lagonly
		solution = lagrangesolve(lgraph, max_iterations=max_iterations_lag, ϵ=ϵ, rel_gap=rel_gap, initialmultipliers=lag_initialmultipliers, lagrangeheuristic=lagrangeheuristic,  maxnoimprove = lag_maxnoimprove, δ = lag_δ, α = lag_α, cpbound = lag_cpbound , user_LB=user_LB,user_UB=user_UB )
		results[:lagrangean_time] += solution.clocktime[end]
		if !lgraph.attributes[:is_infeasible]
			if !haskey(bgraph.attributes,:preprocessed)
				solution = bendersolve(bgraph, cuts=benders_cuts, max_iterations=1, is_nonconvex=is_nonconvex, timelimit=benders_timelimit, LB=benders_LB, ϵ=ϵ, rel_gap=rel_gap, UBupdatefrequency=benders_UBupdatefrequency, verbose=benders_verbose)
			end
		else
			results[:LB] = +Inf 
			results[:UB] = 1e8
			results[:best_avg_x] = []
			results[:best_feasible_x] = []
			return results
		end
		results[:best_avg_x] = lgraph.attributes[:best_avg_x]
		results[:LB] = lgraph.attributes[:LB]
		results[:UB] = lgraph.attributes[:UB]
		results[:best_feasible_x] = lgraph.attributes[:best_feasible_x]
		results[:best_solution_source] = :lgraph
	end


	return results
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
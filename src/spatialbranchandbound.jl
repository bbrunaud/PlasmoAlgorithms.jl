mutable struct BABnode
	lgraph::Plasmo.PlasmoGraph
	bgraph::Plasmo.PlasmoGraph
	LBq
	UBq
	xlbq
	xubq
	best_feasible_x
	best_avg_x
	fathom_type
end

BABnode() = BABnode(lgraph, bgraph, -1e6, 1e6, [], [], [], [], :notfathomed)

function bab_solve(rootnode::BABnode;
	time_limit = 1e4,
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
	lag_cpbound=1e6)
	node_list = []
	active_nodes_indices = []
	fathomed_nodes_indices = []
	solution=crosssolve(rootnode.bgraph, rootnode.lgraph, heuristic=heuristic, max_iterations_lag=max_iterations_lag, max_iterations_benders=max_iterations_benders, benders_cuts=benders_cuts, ϵ=ϵ, rel_gap= rel_gap, benders_UBupdatefrequency=benders_UBupdatefrequency,benders_timelimit=benders_timelimit, benders_verbose=benders_verbose, benders_LB=benders_LB, is_nonconvex=is_nonconvex, lag_α=lag_α, lagrangeheuristic=lagrangeheuristic, lag_initialmultipliers=lag_initialmultipliers, lag_δ=lag_δ, lag_maxnoimprove=lag_maxnoimprove, lag_cpbound=lag_cpbound)
	rootnode.LBq = solution[:LB]
	rootnode.UBq = solution[:UB]
	rootnode.best_feasible_x = deepcopy(solution[:best_feasible_x])
	rootnode.best_avg_x = deepcopy(solution[:best_avg_x])
	LB = solution[:LB]
	UB = solution[:UB]
	solution_time = Dict()
	solution_time[:benders_time] = solution[:benders_time]
	solution_time[:lagrangean_time] = solution[:lagrangean_time]

	xlb_root = Dict()
	xub_root = Dict()
	for varname in keys(rootnode.bgraph.attributes[:roots][1].attributes[:varname_to_var])
		xlb_root[varname] = getlowerbound(rootnode.bgraph.attributes[:roots][1].attributes[:varname_to_var][varname])
		xub_root[varname] = getupperbound(rootnode.bgraph.attributes[:roots][1].attributes[:varname_to_var][varname])
	end

	rootnode.xlbq = xlb_root
	rootnode.xubq = xub_root


	push!(node_list, rootnode)
	push!(active_nodes_indices, 1)
	if calculate_gap(rootnode.LBq, rootnode.UBq) < 0.001
		active_nodes_indices =[]
	end
	N = 3
	num_iter = 1

	while length(active_nodes_indices)>0
		#terminate if time out 
		if solution_time[:benders_time] + solution_time[:lagrangean_time] > time_limit
			break
		end

		#select the node with the lowest lower bound 
		index_selected = 1 
		temp_LB = +Inf
		for index in active_nodes_indices
			if node_list[index].LBq < temp_LB
				temp_LB = node_list[index].LBq 
				index_selected = index
			end
		end
		node_to_branch = node_list[index_selected]
		#remove selected node from active nodes 
		active_nodes_indices = filter(x->x≠ index_selected,active_nodes_indices)
		

		# create two new nodes 
		lchild = copy_node(node_to_branch, UB, benders_cuts=benders_cuts)
		rchild = copy_node(node_to_branch, UB, benders_cuts=benders_cuts)

		#apply branching rule 
		if num_iter %3 == 0
			lchild, rchild = largest_rel_diameter(node_to_branch, lchild, rchild, xlb_root, xub_root)
		else
			lchild, rchild = largest_distance(node_to_branch, lchild, rchild, xlb_root, xub_root, node_to_branch.best_avg_x)
		end

		#set fathom type of node selected 
		node_to_branch.fathom_type = :bybranch

		#solve the two new nodes 
		solve_node(lchild, solution_time, heuristic=heuristic, max_iterations_lag=max_iterations_lag, max_iterations_benders=max_iterations_benders, benders_cuts=benders_cuts, ϵ=ϵ, rel_gap=rel_gap, benders_UBupdatefrequency=benders_UBupdatefrequency,benders_timelimit=benders_timelimit, benders_verbose=benders_verbose, benders_LB=benders_LB, is_nonconvex=is_nonconvex, lag_α=lag_α, lagrangeheuristic=lagrangeheuristic, lag_initialmultipliers=lag_initialmultipliers, lag_δ=lag_δ, lag_maxnoimprove=lag_maxnoimprove, lag_cpbound=lag_cpbound)
		solve_node(rchild, solution_time, heuristic=heuristic, max_iterations_lag=max_iterations_lag, max_iterations_benders=max_iterations_benders, benders_cuts=benders_cuts, ϵ=ϵ, rel_gap=rel_gap, benders_UBupdatefrequency=benders_UBupdatefrequency,benders_timelimit=benders_timelimit, benders_verbose=benders_verbose, benders_LB=benders_LB, is_nonconvex=is_nonconvex, lag_α=lag_α, lagrangeheuristic=lagrangeheuristic, lag_initialmultipliers=lag_initialmultipliers, lag_δ=lag_δ, lag_maxnoimprove=lag_maxnoimprove, lag_cpbound=lag_cpbound)

		#add two new nodes 
		push!(node_list, lchild)
		push!(node_list, rchild)
		push!(active_nodes_indices, length(node_list))
		push!(active_nodes_indices, length(node_list)-1)

		#update bounds 
		LB = +Inf
		for index in active_nodes_indices
			if node_list[index].LBq < LB 
				LB = node_list[index].LBq
			end 
			if node_list[index].UBq < UB 
				UB = node_list[index].UBq
			end
		end

		#check if the two new nodes can be fathomed by optimality 
		if calculate_gap(lchild.LBq, lchild.UBq) < 0.001
			lchild.fathom_type = :byoptimality
			active_nodes_indices = filter!(x->x≠ (length(node_list)-1), active_nodes_indices)
		end 
		if calculate_gap(rchild.LBq, rchild.UBq) <0.001
			rchild.fathom_type = :byoptimality
			active_nodes_indices = filter(x->x≠ length(node_list), active_nodes_indices)
		end


		#fathom nodes by bound 
		for index in active_nodes_indices
			if calculate_gap(node_list[index].LBq, UB) < 0.001
				node_list[index].fathom_type = :bybound
				active_nodes_indices = filter(x->x≠ index, active_nodes_indices)
			end
		end



		num_iter += 1

	end	
	results = Dict()
	results[:node_list] = node_list
	results[:active_nodes_indices] = active_nodes_indices
	results[:solution_time] = solution_time
	results[:LB] = LB 
	results[:UB] = UB 

	return results
	
end

function largest_rel_diameter(node_to_branch::BABnode, lchild::BABnode, rchild::BABnode, xlb_root, xub_root)
	#normalization factor for binary variables
	δ = 0.1

	#select the variable with largest relative diameter 
	temp_diameter = 0
	var_selected = " "
	varname_to_varcat = node_to_branch.bgraph.attributes[:roots][1].attributes[:varname_to_varcat]
	varname_to_var = node_to_branch.bgraph.attributes[:roots][1].attributes[:varname_to_var]
	var_lb = 0
	var_ub = 0 
	for varname in keys(varname_to_varcat)
		lb = getlowerbound(varname_to_var[varname])
		ub = getupperbound(varname_to_var[varname])
		if varname_to_varcat[varname] == :Bin 
			if ((ub-lb) /(xub_root[varname] - xlb_root[varname]) * δ )> temp_diameter
				temp_diameter = ((ub-lb) /(xub_root[varname] - xlb_root[varname])) * δ
				var_selected = varname
				var_lb = lb 
				var_ub = ub 
			end
		else
			if ((ub-lb) /(xub_root[varname] - xlb_root[varname]))  > temp_diameter
				temp_diameter = ((ub-lb) /(xub_root[varname] - xlb_root[varname])) 
				var_selected = varname
				var_lb = lb 
				var_ub = ub 
			end
		end
	end

	#branch on var selected 
	if varname_to_varcat[var_selected] == :Bin
		set_bounds(lchild, var_selected, 0, 0)
		set_bounds(rchild, var_selected, 1, 1)
		lchild.xubq[var_selected] = 0
		rchild.xlbq[var_selected] = 1
	else
		set_bounds(lchild, var_selected, (var_lb+var_ub)/2, var_lb)
		set_bounds(rchild, var_selected, (var_ub), (var_lb+var_ub)/2)
		lchild.xubq[var_selected] = (var_lb+var_ub)/2
		rchild.xlbq[var_selected] = (var_lb+var_ub)/2
	end
	return lchild, rchild 
end

function largest_distance(node_to_branch::BABnode, lchild::BABnode, rchild::BABnode, xlb_root, xub_root, opt_x)
	#normalization factor for binary variables
	δ = 0.1

	#select the variable with largest distance 
	temp_distance = 0
	var_selected = " "
	varname_to_varcat = node_to_branch.bgraph.attributes[:roots][1].attributes[:varname_to_varcat]
	varname_to_var = node_to_branch.bgraph.attributes[:roots][1].attributes[:varname_to_var]
	var_lb = 0
	var_ub = 0 
	for varname in keys(varname_to_varcat)
		lb = getlowerbound(varname_to_var[varname])
		ub = getupperbound(varname_to_var[varname])
		distance = 0.0 
		if opt_x[varname] - lb > ub - opt_x[varname]
			distance = ub - opt_x[varname]
		else
			distance = opt_x[varname] - lb 
		end

		if varname_to_varcat[varname] == :Bin 
			if distance /(xub_root[varname] - xlb_root[varname]) * δ > temp_distance
				temp_distance = distance /(xub_root[varname] - xlb_root[varname]) * δ
				var_selected = varname
				var_lb = lb 
				var_ub = ub 				
			end
		else
			if ( distance /(xub_root[varname] - xlb_root[varname]))  > temp_distance
				temp_distance = distance /(xub_root[varname] - xlb_root[varname])
				var_selected = varname
				var_lb = lb 
				var_ub = ub 				
			end
		end
	end

	#branch on var selected 
	if varname_to_varcat[var_selected] == :Bin
		set_bounds(lchild, var_selected, 0, 0)
		set_bounds(rchild, var_selected, 1, 1)
		lchild.xubq[var_selected] = 0
		rchild.xlbq[var_selected] = 1
	else
		set_bounds(lchild, var_selected, opt_x[var_selected], var_lb)
		set_bounds(rchild, var_selected, (var_ub), opt_x[var_selected])
		lchild.xubq[var_selected] = opt_x[var_selected]
		rchild.xlbq[var_selected] = opt_x[var_selected]
	end
	return lchild, rchild 
end

function set_bounds(node::BABnode, varname, ub, lb)
	for n in values(getnodes(node.bgraph))
		setupperbound(n.attributes[:varname_to_var][varname], ub)
		setlowerbound(n.attributes[:varname_to_var][varname], lb)
		if haskey(n.attributes, :nonconvex_varname_to_var)
			setupperbound(n.attributes[:nonconvex_varname_to_var][varname], ub)
			setlowerbound(n.attributes[:nonconvex_varname_to_var][varname], lb)
		end
	end

	for n in values(getnodes(node.lgraph))
		setupperbound(n.attributes[:varname_to_var][varname], ub)
		setlowerbound(n.attributes[:varname_to_var][varname], lb)
	end
end 

function solve_node(node::BABnode, solution_time; heuristic=:lagfirst, 
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
	lag_cpbound=1e6)
	#check if the best_feasible_x is still feasible in the current node 
	for varname in keys(node.best_feasible_x)
		if node.best_feasible_x[varname] > node.xubq[varname] || node.best_feasible_x[varname] < node.xlbq[varname]
			node.lgraph.attributes[:UB] = +Inf
			node.bgraph.attributes[:UB] = +Inf		
		end
	end


	#change α in lgraph to its initial value 
	node.lgraph.attributes[:α] = [node.lgraph.attributes[:α][1]]
	node.bgraph.attributes[:stalled] = false 
	for n in values(getnodes(node.lgraph))
		n.attributes[:Zsl] = []
		n.attributes[:μ] = []
	end
	solution=crosssolve(node.bgraph, node.lgraph, heuristic=heuristic, max_iterations_lag=max_iterations_lag, max_iterations_benders=max_iterations_benders, benders_cuts=benders_cuts, ϵ=ϵ, rel_gap= rel_gap, benders_UBupdatefrequency=benders_UBupdatefrequency,benders_timelimit=benders_timelimit, benders_verbose=benders_verbose, benders_LB=benders_LB, is_nonconvex=is_nonconvex, lag_α=lag_α, lagrangeheuristic=lagrangeheuristic, lag_initialmultipliers=lag_initialmultipliers, lag_δ=lag_δ, lag_maxnoimprove=lag_maxnoimprove, lag_cpbound=lag_cpbound)
	node.LBq = solution[:LB]
	node.UBq = solution[:UB]
	node.best_feasible_x = deepcopy(solution[:best_feasible_x])
	node.best_avg_x = deepcopy(solution[:best_avg_x])

	solution_time[:benders_time] += solution[:benders_time]
	solution_time[:lagrangean_time] += solution[:lagrangean_time]

end

function copy_node(node::BABnode, UB; benders_cuts=[:LP])
	newnode = deepcopy(node)
	for n in newnode.bgraph.attributes[:roots]
		model = getmodel(n)
		model.internalModel = Void
	end
	newnode.lgraph.attributes[:global_UB] = UB 
	newnode.bgraph.attributes[:global_UB] = UB

	#reset the bounds of CGLP if :LIFT in cuts 
	if :LIFT in benders_cuts
		graph = node.bgraph
    	for index in 1:length(graph.nodes)
      		node = graph.nodes[index]
      		node.attributes[:LP_stalled] = false 
      		node.attributes[:stalled] = false
      		if in_degree(graph, node) != 0 
        		setCGLP(node, graph)
      		end
    	end
	end

	return newnode
end






















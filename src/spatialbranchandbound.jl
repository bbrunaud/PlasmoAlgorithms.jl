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

function bab_solve(rootnode::Plasmo.PlasmoGraph)
	active_nodes = []
	fathomed_nodes = []
	crosssolve(rootnode.bgraph, rootnode.lgraph, max_iterations_lag=3, max_iterations_benders=6, is_nonconvex=true, lag_α = 1.5, lag_δ = 0.3)
	rootnode.LBq = rootnode.bgraph.attributes[:LB]
	LB = rootnode.bgraph.attributes[:LB]
	UB = +Inf
	if rootnode.lgraph.attributes[:UB] < rootnode.bgraph.attributes[:UB]
		UB = rootnode.lgraph.attributes[:UB]
		rootnode.UBq = rootnode.lgraph.attributes[:UB]
		rootnode.best_feasible_x = deepcopy(rootnode.lgraph.attributes[:best_feasible_x])
	else 
		UB = rootnode.bgraph.attributes[:UB]
		rootnode.UBq = rootnode.bgraph.attributes[:UB]
		rootnode.best_feasible_x = deepcopy(rootnode.bgraph.attributes[:best_feasible_x])
	end

	xlb_root = Dict()
	xub_root = Dict()
	for varname in keys(rootnode.bgraph.attributes[:roots][1].attributes[:varname_to_var])
		xlb_root[varname] = getlowerbound(rootnode.bgraph.attributes[:roots][1].attributes[:varname_to_var][varname])
		xub_root[varname] = getupperbound(rootnode.bgraph.attributes[:roots][1].attributes[:varname_to_var][varname])
	end

	rootnode.xlbq = xlb_root
	rootnode.xubq = xub_root




	push!(active_nodes, rootnode)
	# while 
	# solution = crosssolve(node.bgraph, node.lgraph, max_iterations_lag=30, max_iterations_benders=60, is_nonconvex=true, lag_α = 1.5, lag_δ = 0.3)
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

function solve_node(node::BABnode, solution_time)
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
	solution=crosssolve(node.bgraph, node.lgraph, max_iterations_lag=10, max_iterations_benders=60, is_nonconvex=true)
	node.LBq = solution[:LB]
	node.UBq = solution[:UB]
	node.best_feasible_x = deepcopy(solution[:best_feasible_x])
	node.best_avg_x = deepcopy(solution[:best_avg_x])

	solution_time[:benders_time] += solution[:benders_time]
	solution_time[:lagrangean_time] += solution[:lagrangean_time]

end

function copy_node(node::BABnode, UB)
	newnode = deepcopy(node)
	for n in newnode.bgraph.attributes[:roots]
		model = getmodel(n)
		model.internalModel = Void
	end
	newnode.lgraph.attributes[:global_UB] = UB 
	newnode.bgraph.attributes[:global_UB] = UB
	return newnode
end






















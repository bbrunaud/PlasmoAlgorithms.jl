mutable struct BABnode
	lgraph::Plasmo.PlasmoGraph
	bgraph::Plasmo.PlasmoGraph
	LBq
	UBq
	xlbq
	xubq
	best_feasible_x
	fathom_type
end

BABnode() = BABnode(lgraph, bgraph, -1e6, 1e6, [], [], [], :notfathomed)

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
	set_bounds(lchild, var_selected, (var_lb+var_ub)/2, var_lb)
	set_bounds(rchild, var_selected, (var_ub), (var_lb+var_ub)/2)
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

function solve_node(node::BABnode)
	#change α in lgraph to its initial value 
	node.lgraph.attributes[:α] = [node.lgraph.attributes[:α][1]]
	for n in values(getnodes(node.lgraph))
		n.attributes[:Zsl] = []
		n.attributes[:μ] = []
	end
	bendersolve(node.bgraph, max_iterations=6, timelimit=100000, is_nonconvex=true)
	# solution=crosssolve(node.bgraph, node.lgraph, max_iterations_lag=3, max_iterations_benders=3, is_nonconvex=true)
	# node.LBq = solution[:LB]
	# node.UBq = solution[:UB]
	# node.best_feasible_x = deepcopy(solution[:best_feasible_x])

end























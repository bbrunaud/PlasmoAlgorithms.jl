using PlasmoAlgorithms
using Plasmo
using CPLEX
using BARON
using Gurobi
include("input.jl")
include("benderssub.jl")
include("bendersmaster.jl")
include("lagsub.jl")
include("ubsub.jl")


master = generate_bendersmaster()
g = PlasmoGraph()
g.solver = CplexSolver(CPX_PARAM_SCRIND=0)
n1 = add_node(g)
setmodel(n1, master)
println(master)
println(master.colNames)
for s in scenarios
	benderssub = generate_benderssub(psi_f = psi_f, psi_d1=psi_d1, psi_d2=psi_d2, psi_b1=psi_b1, psi_b2=psi_b2,y_up = y_up[:, :, s], DU=DU[:, s], prob=prob[s])
	ubsub = generate_ubsub(psi_f = psi_f, psi_d1=psi_d1, psi_d2=psi_d2, psi_b1=psi_b1, psi_b2=psi_b2, DU=DU[:, s], prob=prob[s])
	n = add_node(g)
	n.attributes[:scenario] = s
	setmodel(n, benderssub)
	add_edge(g, n1, n)
	A = getindex(n1.model,:A)
	@linkconstraint(g, [i in feeds], n1[:A][i] == n[:A][i])
	gamma_intlt = getindex(n1.model, :gamma_intlt)
	@linkconstraint(g, [i in feeds], n1[:gamma_intlt][i] == n[:gamma_intlt][i])
	gamma_pool = getindex(n1.model, :gamma_pool)
	@linkconstraint(g, [i in pools], n1[:gamma_pool][i] == n[:gamma_pool][i])
	S = getindex(n1.model, :S)
	@linkconstraint(g, [i in pools], n1[:S][i] == n[:S][i])
	n.attributes[:ubsub] = ubsub
	println(benderssub)
	println(ubsub)
end

# bendersolve(g, cuts=[:LP], max_iterations=10, timelimit=100000, is_nonconvex=true)

lag_g = PlasmoGraph()
lag_g.solver = BaronSolver(maxtime=5e4, epsr= 1e-3, CplexLibName = "/opt/ibm/ILOG/CPLEX_Studio127/cplex/bin/x86-64_linux/libcplex1270.so")
lagsub1 = generate_lagsub(psi_f = psi_f, psi_d1=psi_d1, psi_d2=psi_d2, psi_b1=psi_b1, psi_b2=psi_b2, DU=DU[:, 1], prob=prob[1])
n1 = add_node(lag_g)
setmodel(n1, lagsub1)
n1.attributes[:prob] =prob[1]
n1.attributes[:scenario] = 1
println(lagsub1)
for s in scenarios
	if s == length(scenarios)
		break
	end
	ss = s + 1
	n = add_node(lag_g)
	n.attributes[:prob] = prob[ss]
	n.attributes[:scenario] = ss 
	lagsub = generate_lagsub(psi_f = psi_f, psi_d1=psi_d1, psi_d2=psi_d2, psi_b1=psi_b1, psi_b2=psi_b2, DU=DU[:, ss], prob=prob[ss])
	setmodel(n, lagsub)
	@linkconstraint(lag_g, [i in feeds], n1[:A][i] == n[:A][i])
	@linkconstraint(lag_g, [i in feeds], n1[:gamma_intlt][i] == n[:gamma_intlt][i])
	@linkconstraint(lag_g, [i in pools], n1[:gamma_pool][i] == n[:gamma_pool][i])
	@linkconstraint(lag_g, [i in pools], n1[:S][i] == n[:S][i])
	println(lagsub)
end
# lagrangesolve(lag_g, max_iterations=5, lagrangeheuristic=nearest_scenario,  maxnoimprove = 1)
# crosssolve(g, lag_g, max_iterations_lag=2, max_iterations_benders=1, is_nonconvex=true)

# all_α = [2 1.5 1 0.5]
# all_δ = [0.8 0.5 0.3]
# # all_α = [2 ]
# # all_δ = [0.8 ]
# gap = []
# for α in all_α
# 	for δ in all_δ
# 		gg = deepcopy(g)
# 		lagg = deepcopy(lag_g)
# 		solution = crosssolve(gg, lagg, max_iterations_lag=30, max_iterations_benders=60, is_nonconvex=true, lag_α = α, lag_δ = δ)
# 		push!(gap, solution.gap)
# 	end
# end
# solution = crosssolve(g, lag_g, max_iterations_lag=30, max_iterations_benders=60, is_nonconvex=true, lag_α = 1.5, lag_δ = 0.3)
# println(gap)

#start debug spatial bab
rootnode = BABnode(lag_g, g, -1e6, 1e6, [], [], [], [], :notfathomed)
node_list = []
active_nodes_indices = []
fathomed_nodes_indices = []
solution=crosssolve(rootnode.bgraph, rootnode.lgraph, benders_cuts=[:GMI], max_iterations_lag=10, max_iterations_benders=60, is_nonconvex=true, lag_α = 1.5, lag_δ = 0.3)
rootnode.LBq = solution[:LB]
rootnode.UBq = solution[:UB]
rootnode.best_feasible_x = deepcopy(solution[:best_feasible_x])
rootnode.best_avg_x = deepcopy(solution[:best_avg_x])
LB = solution[:LB]
UB = solution[:UB]


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
solution_time = Dict()
solution_time[:benders_time] = 0
solution_time[:lagrangean_time] = 0
while length(active_nodes_indices)>0
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
	lchild = copy_node(node_to_branch, UB)
	rchild = copy_node(node_to_branch, UB)

	#apply branching rule 
	if num_iter %3 == 0
		lchild, rchild = largest_rel_diameter(node_to_branch, lchild, rchild, xlb_root, xub_root)
	else
		lchild, rchild = largest_distance(node_to_branch, lchild, rchild, xlb_root, xub_root, node_to_branch.best_avg_x)
	end

	#set fathom type of node selected 
	node_to_branch.fathom_type = :bybranch

	#solve the two new nodes 
	solve_node(lchild, solution_time)
	solve_node(rchild, solution_time)

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
	if length(node_list) > 10
		break
	end 

end






















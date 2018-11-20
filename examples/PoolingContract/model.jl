include("input3.jl")
function generate_model()
	m = Model(solver=BaronSolver(maxtime=5e4, epsr= 1e-3, CplexLibName = "/opt/ibm/ILOG/CPLEX_Studio127/cplex/bin/x86-64_linux/libcplex1270.so"))

	@variable(m, gamma_intlt[i in feeds], Bin)
	@variable(m, gamma_pool[l in pools], Bin)
	@variable(m, SL[l]<=S[l in pools]<=SU[l])
	@variable(m, AL[i]<=A[i in feeds]<=AU[i])
	@variable(m, θ[w in scenarios])
	@variable(m, y[l in pools, j in products, w in scenarios]>=0)
	@variable(m, z[i in feeds, j in products, w in scenarios]>=0)
	@variable(m, q[i in feeds, l in pools, w in scenarios]>=0)
	@variable(m, uf[i in feeds, w in scenarios], Bin)
	@variable(m, ub[i in feeds, w in scenarios], Bin)
	@variable(m, ud[i in feeds, w in scenarios], Bin)
	@variable(m, ub1[i in feeds, w in scenarios], Bin)
	@variable(m, ub2[i in feeds, w in scenarios], Bin)
	@variable(m, ud1[i in feeds, w in scenarios], Bin)
	@variable(m, ud2[i in feeds, w in scenarios], Bin)
	@variable(m, CT[i in feeds, w in scenarios]>=0)
	@variable(m, CTf[i in feeds, w in scenarios]>=0)
	@variable(m, CTb[i in feeds, w in scenarios]>=0)
	@variable(m, CTd[i in feeds, w in scenarios]>=0)
	@variable(m, Bf[i in feeds, w in scenarios]>=0)
	@variable(m, Bd[i in feeds, w in scenarios]>=0)
	@variable(m, Bd1[i in feeds, w in scenarios]>=0)
	@variable(m, Bd2[i in feeds, w in scenarios]>=0)
	@variable(m, Bd11[i in feeds, w in scenarios]>=0)
	@variable(m, Bd12[i in feeds, w in scenarios]>=0)
	@variable(m, Bb[i in feeds, w in scenarios]>=0)
	@variable(m, Bb1[i in feeds, w in scenarios]>=0)
	@variable(m, Bb2[i in feeds, w in scenarios]>=0)
	@variable(m, stage1cost)

	@constraint(m, f1[i in feeds], AL[i]*gamma_intlt[i]<= A[i])
	@constraint(m, f2[i in feeds], A[i] <= AU[i]*gamma_intlt[i])
	@constraint(m, f3[l in pools], SL[l]* gamma_pool[l]<=S[l])
	@constraint(m, f4[l in pools], S[l]<= SU[l]* gamma_pool[l])

	# @constraint(m, e1[i in feeds, w in scenarios], AL[i]*gamma_intlt[i] <= sum(q[i,l,w]*y[l,j,w] for l in pools for j in products if ((i,l) in Tx && (l,j) in Ty)) + sum(z[i,j,w] for j in products if (i,j) in Tz))
	@NLconstraint(m, e2[i in feeds, w in scenarios], sum(q[i,l,w]*y[l,j,w] for l in pools for j in products if ((i,l) in Tx && (l,j) in Ty))  + sum(z[i,j,w] for j in products if (i,j) in Tz)<= A[i])
	@constraint(m, e3[l in pools, w in scenarios], sum(y[l,j,w] for j in products if (l,j) in Ty) <= S[l])
	# @constraint(m, e4[j in products], DL[j] <= sum(y[l,j] for l in pools if (l,j) in Ty) + sum(z[i,j] for i in feeds))
	@constraint(m, e5[j in products, w in scenarios], sum(y[l,j,w] for l in pools if (l,j) in Ty) + sum(z[i,j,w] for i in feeds if (i,j) in Tz)<= DU[j,w])
	@constraint(m, e6[l in pools, w in scenarios], sum(q[i,l,w] for i in feeds if (i,l) in Tx) == 1 )
	
	@NLconstraint(m, e8[j in products, k in qualities, w in scenarios], PL[j,k] * (sum(y[l,j,w] for l in pools if (l,j) in Ty) + sum(z[i,j,w] for i in feeds if (i,j) in Tz)) <= sum(CC[i,k] * z[i,j,w] for i in feeds if (i,j) in Tz) + sum(CC[i,k]*q[i,l,w]*y[l,j,w] for l in pools for i in feeds if ((l,j) in Ty &&(i,l) in Tx)))
	@NLconstraint(m, e9[j in products, k in qualities, w in scenarios], PU[j,k] * (sum(y[l,j,w] for l in pools if (l,j) in Ty) + sum(z[i,j,w] for i in feeds if (i,j) in Tz)) >= sum(CC[i,k] * z[i,j,w] for i in feeds if (i,j) in Tz) + sum(CC[i,k]*q[i,l,w]*y[l,j,w] for l in pools for i in feeds if ((l,j) in Ty &&(i,l) in Tx)))
	# @constraint(m, e10[i in feeds, l in pools; (i,l) in Tx], q[i,l] <=  gamma_pool[l])
	# @constraint(m, e11[i in feeds, l in pools, w in scenarios; (i,l) in Tx], q[i,l,w] <=  gamma_intlt[i])
	# @constraint(m, e12[l in pools, j in products; (l,j) in Ty], y[l,j] <= y_up[l,j] * gamma_pool[l])

	@NLconstraint(m, c1[i in feeds, w in scenarios], sum(q[i,l,w]*y[l,j,w] for l in pools for j in products if ((i,l) in Tx && (l,j) in Ty))  + sum(z[i,j,w] for j in products if (i,j) in Tz)== Bf[i,w]+Bd[i,w]+Bb[i,w])
	@constraint(m, c2[i in feeds, w in scenarios], AL[i]*uf[i,w]<=Bf[i,w])
	@constraint(m, c3[i in feeds, w in scenarios], AU[i]*uf[i,w]>=Bf[i,w])
	@constraint(m, c4[i in feeds, w in scenarios], AL[i]*ud[i,w]<=Bd[i,w])
	@constraint(m, c5[i in feeds, w in scenarios], AU[i]*ud[i,w]>=Bd[i,w])	
	@constraint(m, c6[i in feeds, w in scenarios], AL[i]*ub[i,w]<=Bb[i,w])
	@constraint(m, c7[i in feeds, w in scenarios], AU[i]*ub[i,w]>=Bb[i,w])
	@constraint(m, c8[i in feeds, w in scenarios], uf[i,w]+ub[i,w]+ud[i,w]<= gamma_intlt[i])
	@constraint(m, c9[i in feeds, w in scenarios], CT[i,w]==CTf[i,w]+CTb[i,w]+CTd[i,w])

	#fixed cost 
	@constraint(m, c10[i in feeds, w in scenarios], CTf[i,w]==psi_f[i,w]*Bf[i,w])

	#dicount after a certain amount
	@constraint(m, c11[i in feeds, w in scenarios], CTd[i,w]==psi_d1[i,w]*Bd1[i,w]+psi_d2[i,w]*Bd2[i,w])
	@constraint(m, c12[i in feeds, w in scenarios], Bd[i,w]== Bd1[i,w]+Bd2[i,w])
	@constraint(m, c13[i in feeds, w in scenarios], Bd1[i,w]== Bd11[i,w]+Bd12[i,w])
	@constraint(m, c14[i in feeds, w in scenarios], Bd11[i,w]<= sigma_d[i]*ud1[i,w])
	@constraint(m, c15[i in feeds, w in scenarios], Bd12[i,w]==sigma_d[i]*ud2[i,w])
	@constraint(m, c16[i in feeds, w in scenarios], Bd2[i,w]<= AU[i]*ud2[i,w])

	#bulk discount contract 
	@constraint(m, c17[i in feeds, w in scenarios], CTb[i,w]== psi_b1[i,w]*Bb1[i,w]+psi_b2[i,w]*Bb2[i,w])
	@constraint(m, c18[i in feeds, w in scenarios], Bb[i,w]== Bb1[i,w]+Bb2[i,w])
	@constraint(m, c19[i in feeds, w in scenarios], Bb1[i,w]<= sigma_b[i]*ub1[i,w])
	@constraint(m, c20[i in feeds, w in scenarios], Bb2[i,w]>= sigma_b[i]*ub2[i,w])
	@constraint(m, c21[i in feeds, w in scenarios], Bb2[i,w]<= AU[i]*ub2[i,w])
	@constraint(m, c22[i in feeds, w in scenarios], ub1[i,w]+ub2[i,w]==ub[i,w])

	# @constraint(m, cobj[w in scenarios], θ[w] == prob[w]*(sum(CT[i,w] for i in feeds) -sum(d[j,w]*(sum(y[l,j,w] for l in pools if (l,j) in Ty) + sum(z[i,j,w] for i in feeds if (i,j) in Tz)) for j in products) ))
	@constraint(m, cobj[w in scenarios], θ[w] == (sum(CT[i,w] for i in feeds) -sum(d[j,w]*(sum(y[l,j,w] for l in pools if (l,j) in Ty) + sum(z[i,j,w] for i in feeds if (i,j) in Tz)) for j in products) ))

	@constraint(m, stage1cost == sum(c_fixed_inlt[i] * gamma_intlt[i]+c_variable_inlt[i]*A[i] for i in feeds) + sum(c_fixed_pool[l] * gamma_pool[l] + c_variable_pool[l]*S[l] for l in pools) )
	@objective(m, Min,  sum(c_fixed_inlt[i] * gamma_intlt[i]+c_variable_inlt[i]*A[i] for i in feeds) + sum(c_fixed_pool[l] * gamma_pool[l] + c_variable_pool[l]*S[l] for l in pools) + sum(prob[w]*(sum(CT[i,w] for i in feeds) -sum(d[j,w]*(sum(y[l,j,w] for l in pools if (l,j) in Ty) + sum(z[i,j,w] for i in feeds if (i,j) in Tz)) for j in products) ) for w in scenarios))

	return m
end
# r1 = 1.3
# r2 = 1.5
# r3 =  1.5
# d = base_d * r1
# c_variable_pool = base_c_variable_pool * r2 
# c_variable_inlt = base_c_variable_inlt * r2 
# psi_f = zeros(length(feeds), length(scenarios))
# psi_d1  = zeros(length(feeds), length(scenarios))
# psi_d2 = zeros(length(feeds), length(scenarios))
# psi_b1 = zeros(length(feeds), length(scenarios))
# psi_b2 = zeros(length(feeds), length(scenarios))

# ratio =[0.7 1.0 1.3]
# for w in 1:length(scenarios)
# 	psi_f[:, w] = base_psi_f * ratio[w] * r3
# 	psi_d2[:,w] = base_psi_d2 * ratio[w] * r3
# 	psi_d1[:,w] = base_psi_d1 * ratio[w] * r3
# 	psi_b1[:,w] = base_psi_b1 * ratio[w] * r3
# 	psi_b2[:,w] = base_psi_b2 * ratio[w] * r3 
# end
model = generate_model()
solve(model)
println(getobjectivevalue(model))
println(getvalue(getindex(model, :A)))
println(getvalue(getindex(model, :S)))
println(getvalue(getindex(model, :θ)))
println(getvalue(getindex(model, :stage1cost)))
# println(getvalue(getindex(model, :y)))
# println(getvalue(getindex(model, :z)))
# y_value = getvalue(getindex(model, :y))
# q_value = getvalue(getindex(model, :q))
# x_value = zeros(length(feeds), length(pools), length(scenarios))
# for w in scenarios
# 	for i in feeds
# 		for l in pools
# 			x_value[i,l,w] = sum(y_value[l,j,w]*q_value[i,l,w] for j in products)
# 		end
# 	end
# end
# println("xvalue")
# for w in scenarios
# 	println("scenarios ", w)
# 	for i in feeds
# 		println("feed ", i)
# 		println(x_value[i,:,w])
# 	end
# end

# println(getvalue(getindex(model, :q)))
# println(getvalue(getindex(model, :gamma_intlt)))
# println(getvalue(getindex(model, :gamma_pool)))
# println(getvalue(getindex(model, :uf)))
# println(getvalue(getindex(model, :ub)))
# println(getvalue(getindex(model, :ud)))
# println(getvalue(getindex(model, :Bf)))
# println(getvalue(getindex(model, :Bd)))
# println(getvalue(getindex(model, :Bd1)))
# println(getvalue(getindex(model, :Bd2)))
# println(getvalue(getindex(model, :Bd11)))
# println(getvalue(getindex(model, :Bd12)))
# println(getvalue(getindex(model, :Bb)))
# println(getvalue(getindex(model, :Bb1)))
# println(getvalue(getindex(model, :Bb2)))


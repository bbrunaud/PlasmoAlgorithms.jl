function generate_benderssub(;psi_f=base_psi_f, psi_d1=base_psi_d1, psi_d2=base_psi_d2, psi_b1=base_psi_b1, psi_b2=base_psi_b2, DU=base_DU, prob=0)
	m = Model(solver=CplexSolver())

	@variable(m, gamma_intlt[i in feeds], Bin)
	@variable(m, gamma_pool[l in pools], Bin)
	@variable(m, S[l in pools]>=0)
	@variable(m, A[i in feeds]>=0)
	@variable(m, y[l in pools, j in products]>=0)
	@variable(m, z[i in feeds, j in products]>=0)
	@variable(m, 0<=q[i in feeds, l in pools]<=1)
	@variable(m, uf[i in feeds], Bin)
	@variable(m, ub[i in feeds], Bin)
	@variable(m, ud[i in feeds], Bin)
	@variable(m, ub1[i in feeds], Bin)
	@variable(m, ub2[i in feeds], Bin)
	@variable(m, ud1[i in feeds], Bin)
	@variable(m, ud2[i in feeds], Bin)
	@variable(m, CT[i in feeds]>=0)
	@variable(m, CTf[i in feeds]>=0)
	@variable(m, CTb[i in feeds]>=0)
	@variable(m, CTd[i in feeds]>=0)
	@variable(m, Bf[i in feeds]>=0)
	@variable(m, Bd[i in feeds]>=0)
	@variable(m, Bd1[i in feeds]>=0)
	@variable(m, Bd2[i in feeds]>=0)
	@variable(m, Bd11[i in feeds]>=0)
	@variable(m, Bd12[i in feeds]>=0)
	@variable(m, Bb[i in feeds]>=0)
	@variable(m, Bb1[i in feeds]>=0)
	@variable(m, Bb2[i in feeds]>=0)

	#q * y 
	@variable(m, qy[i in feeds, l in pools, j in products]>=0)

	#discretize q 
	@variable(m, dot_q[i in feeds, l in pools, t in intervals]>=0)
	@variable(m, dot_y[l in pools, j in pools, i in feeds, t in intervals]>=0)
	@variable(m, dot_qy[i in feeds, l in pools, j in products, t in intervals]>=0)
	@variable(m, delta[i in feeds, l in pools, t in intervals], Bin)

	# @constraint(m, e1[i in feeds], AL[i]*gamma_intlt[i] <= sum(q[i,l]*y[l,j] for l in pools for j in products if ((i,l) in Tx && (l,j) in Ty)) + sum(z[i,j] for j in products if (i,j) in Tz))
	@constraint(m, e2[i in feeds], sum(qy[i,l,j] for l in pools for j in products if ((i,l) in Tx && (l,j) in Ty))  + sum(z[i,j] for j in products if (i,j) in Tz)<= A[i])
	@constraint(m, e3[l in pools], sum(y[l,j] for j in products if (l,j) in Ty) <= S[l])
	# @constraint(m, e4[j in products], DL[j] <= sum(y[l,j] for l in pools if (l,j) in Ty) + sum(z[i,j] for i in feeds))
	@constraint(m, e5[j in products], sum(y[l,j] for l in pools if (l,j) in Ty) + sum(z[i,j] for i in feeds if (i,j) in Tz)<= DU[j])
	@constraint(m, e6[l in pools], sum(q[i,l] for i in feeds if (i,l) in Tx) == 1)
	
	@constraint(m, e8[j in products, k in qualities], PL[j,k] * (sum(y[l,j] for l in pools if (l,j) in Ty) + sum(z[i,j] for i in feeds if (i,j) in Tz)) <= sum(CC[i,k] * z[i,j] for i in feeds if (i,j) in Tz) + sum(CC[i,k]*qy[i,l,j] for l in pools for i in feeds if ((l,j) in Ty &&(i,l) in Tx)))
	@constraint(m, e9[j in products, k in qualities], PU[j,k] * (sum(y[l,j] for l in pools if (l,j) in Ty) + sum(z[i,j] for i in feeds if (i,j) in Tz)) >= sum(CC[i,k] * z[i,j] for i in feeds if (i,j) in Tz) + sum(CC[i,k]*qy[i,l,j] for l in pools for i in feeds if ((l,j) in Ty &&(i,l) in Tx)))
	# @constraint(m, e10[i in feeds, l in pools; (i,l) in Tx], q[i,l] <=  gamma_pool[l])
	# @constraint(m, e11[i in feeds, l in pools; (i,l) in Tx], q[i,l] <=  gamma_intlt[i])
	# @constraint(m, e12[l in pools, j in products; (l,j) in Ty], y[l,j] <= y_up[l,j] * gamma_pool[l])

# 	#the piecewise McCormick ====================
	@constraint(m, p1[i in feeds ,l in pools, j in products], qy[i,l,j] == sum(dot_qy[i,l,j,t] for t in intervals))
	@constraint(m, p2[l in pools, j in products, i in feeds], y[l,j] == sum(dot_y[l,j,i,t] for t in intervals))
	@constraint(m, p3[i in feeds, l in pools], q[i,l] == sum(dot_q[i,l,t] for t in intervals))
	@constraint(m, p4[i in feeds, l in pools], sum(delta[i,l,t] for t in intervals) == 1)
	@constraint(m, Mc1[i in feeds, l in pools, j in products, t in intervals], dot_qy[i,l,j,t]>= q_points[t+1] * dot_y[l,j,i,t] + dot_q[i,l,t] * y_up[l,j] - q_points[t+1] * y_up[l,j] * delta[i,l,t] )
	@constraint(m, Mc2[i in feeds, l in pools, j in products, t in intervals], dot_qy[i,l,j,t]>= q_points[t] * dot_y[l,j,i,t]  )
	@constraint(m, Mc3[i in feeds, l in pools, j in products, t in intervals], dot_qy[i,l,j,t]<= q_points[t+1] * dot_y[l,j,i,t]  )
	@constraint(m, Mc4[i in feeds, l in pools, j in products, t in intervals], dot_qy[i,l,j,t]<= q_points[t] * dot_y[l,j,i,t] + dot_q[i,l,t] * y_up[l,j] - q_points[t] * y_up[l,j] * delta[i,l,t] )
	@constraint(m, b1[i in feeds, l in pools, t in intervals], dot_q[i,l,t]<= q_points[t+1]*delta[i,l,t])
	@constraint(m, b2[i in feeds, l in pools, t in intervals], dot_q[i,l,t]>= q_points[t]*delta[i,l,t])
	@constraint(m, b3[l in pools, j in products, i in feeds, t in intervals], dot_y[l,j,i,t] <= y_up[l,j] * delta[i,l,t])
	# #end ======================

	@constraint(m, c1[i in feeds], sum(qy[i,l,j] for l in pools for j in products if ((i,l) in Tx && (l,j) in Ty))  + sum(z[i,j] for j in products if (i,j) in Tz)== Bf[i]+Bd[i]+Bb[i])
	@constraint(m, c2[i in feeds], AL[i]*uf[i]<=Bf[i])
	@constraint(m, c3[i in feeds], AU[i]*uf[i]>=Bf[i])
	@constraint(m, c4[i in feeds], AL[i]*ud[i]<=Bd[i])
	@constraint(m, c5[i in feeds], AU[i]*ud[i]>=Bd[i])	
	@constraint(m, c6[i in feeds], AL[i]*ub[i]<=Bb[i])
	@constraint(m, c7[i in feeds], AU[i]*ub[i]>=Bb[i])
	@constraint(m, c8[i in feeds], uf[i]+ub[i]+ud[i]<= gamma_intlt[i])
	@constraint(m, c9[i in feeds], CT[i]==CTf[i]+CTb[i]+CTd[i])

	#fixed cost 
	@constraint(m, c10[i in feeds], CTf[i]==psi_f[i]*Bf[i])

	#dicount after a certain amount
	@constraint(m, c11[i in feeds], CTd[i]==psi_d1[i]*Bd1[i]+psi_d2[i]*Bd2[i])
	@constraint(m, c12[i in feeds], Bd[i]== Bd1[i]+Bd2[i])
	@constraint(m, c13[i in feeds], Bd1[i]== Bd11[i]+Bd12[i])
	@constraint(m, c14[i in feeds], Bd11[i]<= sigma_d[i]*ud1[i])
	@constraint(m, c15[i in feeds], Bd12[i]==sigma_d[i]*ud2[i])
	@constraint(m, c16[i in feeds], Bd2[i]<= AU[i]*ud2[i])

	#bulk discount contract 
	@constraint(m, c17[i in feeds], CTb[i]== psi_b1[i]*Bb1[i]+psi_b2[i]*Bb2[i])
	@constraint(m, c18[i in feeds], Bb[i]== Bb1[i]+Bb2[i])
	@constraint(m, c19[i in feeds], Bb1[i]<= sigma_b[i]*ub1[i])
	@constraint(m, c20[i in feeds], Bb2[i]>= sigma_b[i]*ub2[i])
	@constraint(m, c21[i in feeds], Bb2[i]<= AU[i]*ub2[i])
	@constraint(m, c22[i in feeds], ub1[i]+ub2[i]==ub[i])

	@objective(m, Min,   prob*(sum(CT[i] for i in feeds) -sum(d[j]*(sum(y[l,j] for l in pools if (l,j) in Ty) + sum(z[i,j] for i in feeds if (i,j) in Tz)) for j in products) ) )

	return m
end
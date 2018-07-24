include("model.jl")

ratio1 = 1:20
ratio1 = ratio1/10
good_capacity_ratio = []
good_scenario_ratio = []
ratio2= [0.5 1 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0]
ratio3 = [0.5 1 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0]
for r1 in ratio1
	d = base_d * r1 
	for r2 in ratio2
		c_variable_pool = base_c_variable_pool * r2 
		c_variable_inlt = base_c_variable_inlt * r2 
		for r3 in ratio3
			psi_f = zeros(length(feeds), length(scenarios))
			psi_d1  = zeros(length(feeds), length(scenarios))
			psi_d2 = zeros(length(feeds), length(scenarios))
			psi_b1 = zeros(length(feeds), length(scenarios))
			psi_b2 = zeros(length(feeds), length(scenarios))

			ratio =[0.7 1.0 1.3]
			for w in 1:length(scenarios)
				psi_f[:, w] = base_psi_f * ratio[w] * r3
				psi_d2[:,w] = base_psi_d2 * ratio[w] * r3
				psi_d1[:,w] = base_psi_d1 * ratio[w] * r3
				psi_b1[:,w] = base_psi_b1 * ratio[w] * r3
				psi_b2[:,w] = base_psi_b2 * ratio[w] * r3 
			end
			println(psi_f)
			m = generate_model()
			solve(m)
			println("===========ratio==================")
			println(r1)
			println(r2)
			println(r3)
			println(getobjectivevalue(m))
			println(getvalue(getindex(m, :A)))
			println(getvalue(getindex(m, :S)))
			# println(getvalue(getindex(m, :y)))
			# println(getvalue(getindex(m, :z)))
			# println(getvalue(getindex(m, :q)))
			# println(getvalue(getindex(m, :gamma_intlt)))
			# println(getvalue(getindex(m, :gamma_pool)))
			# println(getvalue(getindex(m, :uf)))
			# println(getvalue(getindex(m, :ub)))
			# println(getvalue(getindex(m, :ud)))
			println(getvalue(getindex(m, :Bf)))
			println(getvalue(getindex(m, :Bd)))
			println(getvalue(getindex(m, :Bb)))	
			is_capacity = false 
			is_scenario = false
			A = getvalue(getindex(m, :A))
			if sum(A) > 10 && sum(A) < sum(DU) - 10
				is_capacity = true 
				push!(good_capacity_ratio, [r1,r2,r3])
			end
			uf = getvalue(getindex(m, :uf))
			ub = getvalue(getindex(m, :ub))
			ud = getvalue(getindex(m, :ud))
			for i in feeds
				for w in 1:(length(scenarios)-1)
					if uf[i,1] != uf[i, w+1]
						is_scenario = true 
					end
					if ub[i,1] != ub[i, w+1]
						is_scenario = true 
					end
					if ud[i,1] != ud[i, w+1]
						is_scenario = true 
					end	
				end
			end
			if is_scenario
				push!(good_scenario_ratio, [r1,r2,r3])
			end
		end
	end
end

println(good_capacity_ratio)
println(good_scenario_ratio)

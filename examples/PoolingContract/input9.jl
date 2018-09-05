using JuMP
using BARON
#sets
feeds = 1:5
products  = 1:3
pools = 1:4
qualities = 1:2
intervals = 1:1
points = 1:2

q_points = zeros(length(points))
for i in intervals
	q_points[i+1] = (1.0 * i ) / length(intervals)
end

Tx = ((1,1),(1,2),(1,3),(1,4),(2,1),(2,2),(2,3),(2,4),(3,1),(3,2),(3,3),(3,4),(4,1),(4,2),(4,3),(4,4),(5,1),(5,2),(5,3),(5,4))
Ty = ((1,1),(1,2),(1,3),(2,1),(2,2),(2,3),(3,1),(3,2),(3,3),(4,1),(4,2),(4,3))
Tz = ()
# Tz = ((1,1),(1,2),(1,3),(2,1),(2,2),(2,3),(3,1),(3,2),(3,3),(4,1),(4,2),(4,3),(5,1),(5,2),(5,3))
AU = [300 250 250 200 300] 	## Maximum Available Flow
AL = zeros(length(feeds)) 	## Minimum Available Flow
SL = zeros(length(pools))
SU = [400 400 350 500] # # Pool Size
base_d = [5.7  6.2  6.8]	## Product Unit Price
	

base_DU = [229 173 284] 	## Maximum Product Demand	
c_fixed_pool = [310 470 380 510]## Cost Parameter
base_c_variable_pool = [1.1 0.9 1.05 0.8]
c_variable_pool = [1.1 0.9 1.05 0.8]
c_fixed_inlt = [260 70 150 190 110] ## Cost Parameter
c_variable_inlt = [0.5 0.8 0.6 0.55 0.7] 
base_c_variable_inlt = [0.5 0.8 0.6 0.55 0.7] 
CC = [ 		## Feed Concentration
   0.13  0.87
   0.89  0.11 
   0.69  0.31
   0.28  0.72
   0.35  0.65]
      
PU	 =[	## Maximum Allowable Product Concentration
  0.56  0.44
  0.30  0.70
  0.41  0.59]

PL =[	        ## Minimum Allowable Product Concentration
0.56  0.44
0.30  0.70
0.41  0.59]

c_x =[	## Cost Parameter
6.2     9.4     7.6     10.2
1.67    2.53    2.05    2.75
3.58    5.42    4.39    5.89
4.53    6.87    5.55    7.45
2.62    3.98    3.22    4.32]




scenarios=1:9
base_prob = [0.3 0.4 0.3]
prob = zeros(length(scenarios))

sigma_d = AU /2
sigma_b =  AU /3 * 2
base_psi_f = ones(length(feeds))*0.5
base_psi_d1 = ones(length(feeds))*0.55
base_psi_d2 = ones(length(feeds))*0.4
base_psi_b1 = ones(length(feeds))*0.55
base_psi_b2 = ones(length(feeds))*0.48	## Feed Unit Price

psi_f = zeros(length(feeds), length(scenarios))
psi_d1  = zeros(length(feeds), length(scenarios))
psi_d2 = zeros(length(feeds), length(scenarios))
psi_b1 = zeros(length(feeds), length(scenarios))
psi_b2 = zeros(length(feeds), length(scenarios))
d = zeros(length(products), length(scenarios))
DU = zeros(length(products), length(scenarios))
ratio =[0.7 1.0 1.3]
for w1 in 1:3
	for w2 in 1:3
		w = (w2-1) * 3 + w1 
		prob[w] = base_prob[w1] * base_prob[w2]
		DU[:,w] = base_DU * ratio[w1]
		d[:,w] = base_d 
		psi_f[:, w] = base_psi_f * ratio[w2]
		psi_d2[:,w] = base_psi_d2 * ratio[w2]
		psi_d1[:,w] = base_psi_d1 * ratio[w2]
		psi_b1[:,w] = base_psi_b1 * ratio[w2]
		psi_b2[:,w] = base_psi_b2 * ratio[w2]
	end
end


#bounds
x_up = zeros(length(feeds), length(pools), length(scenarios))
for i in feeds
	for l in pools
		for w in scenarios
			x_up[i,l,w] = min(AU[i], SU[l], sum(DU[j,w] for j in products))
		end
	end
end

y_up = zeros(length(pools), length(products), length(scenarios))
for l in pools
	for j in products
		for w in scenarios
			y_up[l,j,w] = min(SU[l], DU[j,w], sum(AU[i] for i in feeds))
		end
	end
end

z_up = zeros(length(feeds), length(products), length(scenarios))
for i in feeds
	for j in products
		for w in scenarios
			z_up[i,j,w] = min(AU[i], DU[j,w])
		end
	end
end

p_lo = zeros(length(pools), length(qualities))
p_up = zeros(length(pools),length(qualities))
for l in pools
	for k in qualities
		p_lo[l,k] = min(CC[:,k]...)
		p_up[l,k] = max(CC[:,k]...)
	end
end








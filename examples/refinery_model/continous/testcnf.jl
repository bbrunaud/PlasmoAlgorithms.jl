using JuMP
using CPLEX
using PlasmoAlgorithms
#sets
feeds = 1:5
products  = 1:3
pools = 1:4
qualities = 1:1
Tx = ((1,1),(1,2),(1,3),(1,4),(2,1),(2,2),(2,3),(2,4),(3,1),(3,2),(3,3),(3,4),(4,1),(4,2),(4,3),(4,4),(5,1),(5,2),(5,3),(5,4))
Ty = ((1,1),(1,2),(1,3),(2,1),(2,2),(2,3),(3,1),(3,2),(3,3),(4,1),(4,2),(4,3))
Tz = ()

d = [10.5 11.2 12.5]    ## Product Unit Price
c = [0  0  0 0  0   ]   ## Feed Unit Price
AU = ones(length(feeds)) * 1e6  ## Maximum Available Flow
AL = zeros(length(feeds))   ## Minimum Available Flow
SS = ones(length(pools)) * 1e6   # # Pool Size
DU = [229 173 284]  ## Maximum Product Demand
DL = [229 173 284]  ## Minimum Product Demand
c_gamma_pool = [310 470 380 510]## Cost Parameter
c_gamma_inlt = [260 70 150 190 110] ## Cost Parameter
CC = [      ## Feed Concentration
   0.13  0.87
   0.89  0.11 
   0.69  0.31
   0.28  0.72
   0.35  0.65]
      
PU   =[ ## Maximum Allowable Product Concentration
  0.56  0.44
  0.30  0.70
  0.41  0.59]

PL =[           ## Minimum Allowable Product Concentration
0.56  0.44
0.30  0.70
0.41  0.59]

c_x =[  ## Cost Parameter
6.2     9.4     7.6     10.2
1.67    2.53    2.05    2.75
3.58    5.42    4.39    5.89
4.53    6.87    5.55    7.45
2.62    3.98    3.22    4.32]

#bounds
x_up = zeros(length(feeds), length(pools))
for i in feeds
    for l in pools
        x_up[i,l] = min(AU[i], SS[l], sum(DU[j] for j in products))
    end
end

y_up = zeros(length(pools), length(products))
for l in pools
    for j in products
        y_up[l,j] = min(SS[l], DU[j], sum(AU[i] for i in feeds))
    end
end

z_up = zeros(length(feeds), length(products))
for i in feeds
    for j in products
        z_up[i,j] = min(AU[i], DU[j])
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

function generate_model()
    model = Model(solver=CplexSolver())
    @variable(model, 0<=y[l in pools, j in products]<=y_up[l,j])
    @variable(model, z[i in feeds, j in products]>=0)
    @variable(model, 0<=q[i in feeds, l in pools]<=1)
    @variable(model, gamma_intlt[i in feeds], Bin)
    @variable(model, gamma_pool[l in pools], Bin)
    @variable(model, qy[i in feeds, l in pools, j in products]>=0)

    @constraint(model, e1[i in feeds], AL[i]*gamma_intlt[i] <= sum(qy[i,l,j] for l in pools for j in products if ((i,l) in Tx && (l,j) in Ty)) + sum(z[i,j] for j in products if (i,j) in Tz))
    @constraint(model, e2[i in feeds], sum(qy[i,l,j] for l in pools for j in products if ((i,l) in Tx && (l,j) in Ty))  + sum(z[i,j] for j in products if (i,j) in Tz)<= AU[i] * gamma_intlt[i])
    @constraint(model, e3[l in pools], sum(y[l,j] for j in products if (l,j) in Ty) <= SS[l] * gamma_pool[l])
    @constraint(model, e4[j in products], DL[j] <= sum(y[l,j] for l in pools if (l,j) in Ty) + sum(z[i,j] for i in feeds))
    @constraint(model, e5[j in products], sum(y[l,j] for l in pools if (l,j) in Ty) + sum(z[i,j] for i in feeds if (i,j) in Tz)<= DU[j])
    @constraint(model, e6[l in pools], sum(q[i,l] for i in feeds if (i,l) in Tx) == gamma_pool[l] )
    
    @constraint(model, e8[j in products, k in qualities], PL[j,k] * (sum(y[l,j] for l in pools if (l,j) in Ty) + sum(z[i,j] for i in feeds if (i,j) in Tz)) <= sum(CC[i,k] * z[i,j] for i in feeds if (i,j) in Tz) + sum(CC[i,k]*qy[i,l,j] for l in pools for i in feeds if ((l,j) in Ty &&(i,l) in Tx)))
    @constraint(model, e9[j in products, k in qualities], PU[j,k] * (sum(y[l,j] for l in pools if (l,j) in Ty) + sum(z[i,j] for i in feeds if (i,j) in Tz)) >= sum(CC[i,k] * z[i,j] for i in feeds if (i,j) in Tz) + sum(CC[i,k]*qy[i,l,j] for l in pools for i in feeds if ((l,j) in Ty &&(i,l) in Tx)))
    @constraint(model, e10[i in feeds, l in pools; (i,l) in Tx], q[i,l] <=  gamma_pool[l])
    @constraint(model, e11[i in feeds, l in pools; (i,l) in Tx], q[i,l] <=  gamma_intlt[i])
    @constraint(model, e12[l in pools, j in products; (l,j) in Ty], y[l,j] <= y_up[l,j] * gamma_pool[l])
    @constraint(model, e13[l in pools, j in products], sum(qy[i,l,j] for i  in feeds if((i,l) in Tx)) == y[l,j])
    @objective(model, Min, sum((c_x[i, l] + c[i])*qy[i,l,j] for i in feeds for j in products for l in pools if ( (i,l) in Tx && (l,j) in Ty)) - sum(d[j] * y[l,j] for l in pools for j in products if (l,j) in Ty) + sum((-d[j]+c[i])*z[i,j] for i in feeds for j in products if (i,j) in Tz) + sum(c_gamma_inlt[i] * gamma_intlt[i] for i in feeds) + sum(c_gamma_pool[l] * gamma_pool[l] for l in pools))


    A_lp,b_lp,m_lp, n = getModelInfo(model)


    npoints = 5
    for i in feeds
        for l in pools
            add_disjunction(model, A_lp,b_lp,m_lp, n, q[i,l], y[l, :], qy[i,l,:], npoints)
        end
    end
    return model

end 


model = generate_model()
solve(model)
# println(model)






























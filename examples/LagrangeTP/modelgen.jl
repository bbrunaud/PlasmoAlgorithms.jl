using JuMP
using Gurobi

include("ex1inputs.jl")

    oproducts = 1:3
    #m
    markets = 1:3
    #s
    sites = 1:3
    #t
    otime = 1:3
    #nt

function spmodel(products=oproducts,time=otime)

    m1 = Model(solver=GurobiSolver(OutputFlag=0))

###Positive Variables###
#Sales
@variable(m1, sl[ m in markets,i in products, t in time] >= 0)
#Inventory to site
@variable(m1, vi[s in sites, i in products, t in time] >= 0)
@variable(m1, vf[s in sites, i in products, t in time] >= 0)
#Production
@variable(m1, x[s in sites, i in products, t in time] >= 0)
#Shipment received
@variable(m1, y[m in markets, s in sites, i in products, t in time] >= 0)
#Time spent producing i in s during t
@variable(m1, tt[s in sites, i in products, t in time] >= 0)

#Decompostion Variables
@variable(m1, hi[s in sites, i in products, t in time] >= 0, upperbound=HT)
#Dual Variable
@variable(m1, hf[s in sites, i in products, t in time] >= 0, upperbound=HT)

###Binary Variables###
#1 if site s produces product i in time period t
@variable(m1, setup[s in sites, i in products, t in time], Bin)

###Equations###
#Site Mass Balance
@constraint(m1, eq1b[s in sites, i in products, t in time], x[s,i,t]+ vi[s,i,t] == sum(y[m,s,i,t] for m in markets)+vf[s,i,t])
if 1 in time
  @constraint(m1, eq1a[s in sites, i in products], vi[s,i,1] == 0)
end

#a descomponer
@constraint(m1, eq3a[s in sites, i in products, t in time], tt[s,i,t]+setup[s,i,t]*sut[s,i] <= hi[s,i,t])
@constraint(m1, eq3b[s in sites,i in products, t in time], hf[s,i,t]==hi[s,i,t]-(tt[s,i,t]+setup[s,i,t]*sut[s,i]))
#dual
if 1 in products
  @constraint(m1, eq3c[s in sites, t in time], hi[s,1,t]==HT)
end

#Market mass balances
@constraint(m1, eq4[m in markets, i in products, t in time], sl[m,i,t] == sum(y[m,s,i,t] for s in sites))
@constraint(m1, eq5[m in markets, i in products, t in time], sl[m,i,t] <= fcast[m,i,t])

#Shipment mass balance
@constraint(m1, eq7[s in sites, i in products, t in time], x[s,i,t] == tt[s,i,t]*prate[s,i])

@constraint(m1, newconstraint[s in sites, i in products, t in time], tt[s,i,t] <= setup[s,i,t]*HT)
#Objective Function
@expression(m1, sales, sum(β[m,i]*sl[m,i,t] for m in markets, i in products, t in time))
@expression(m1, opercost, sum(α[s,i]*x[s,i,t] for s in sites, i in products, t in time))
@expression(m1, invcost, sum(δ[s,i]*vf[s,i,t] for s in sites,  i in products, t in time))
@expression(m1, shipmentcost, sum(γ[m,s,i]*y[m,s,i,t] for m in markets, s in sites, i in products, t in time))
@expression(m1, septupcosts, sum(setupcost[s,i]*setup[s,i,t] for s in sites, i in products, t in time))

@objective(m1, Max, sales - opercost - invcost - shipmentcost - septupcosts)

###Bounds###
for s in sites, i in products, t in time
    setupperbound(vf[s,i,t], 4*HT*prate[s,i])
    setupperbound(vi[s,i,t], 4*HT*prate[s,i])
end

for m in markets, s in sites, i in products, t in time
    setupperbound(y[m,s,i,t], fcast[m,i,t])
end

# Linkings
@constraint(m1, linkp[s in sites, i in products, t in time; i < length(oproducts) && i+1 in products], hf[s,i,t] == hi[s,i+1,t])
@constraint(m1, linkt[s in sites, i in products, t in time; t < length(otime) && t+1 in time], vf[s,i,t] == vi[s,i,t+1])


return m1
end

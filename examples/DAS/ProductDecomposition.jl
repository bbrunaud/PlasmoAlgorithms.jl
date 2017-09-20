using JuMP
using Gurobi

include("dasset.jl")


cu = ["CUS$i" for i in 1:10]
pr = ["PROD$i" for i in 1:10]
tp = 1:12

grb = GurobiSolver(MIPGap=0.01, Threads=Sys.CPU_CORES)
function model(customers=cu,products=pr,periods=tp;solver=grb)
  m = Model(solver=solver)

  if typeof(products) == String
    products = [products]
  end

  RT = Route
  Dem = Demand
  numProducts = length(products)

  origins = vcat(plants,warehouses)
  destinations = vcat(warehouses,customers)

  M = 1.2*sum(Dem[k,p,t] for k in customers, p in products, t in periods)

  @variables m begin
    x[i in origins, j in destinations, p in products, t in periods; RT[i,j,p]>0] >= 0
    si[j in warehouses, p in products, t in periods] >= 0
    sf[j in warehouses, p in products, t in periods] >= 0
    y[j in warehouses, t in periods, p in products], Bin
    capi[i in plants, p in products, t in periods] >= 0
    capf[i in plants, p in products, t in periods] >= 0
    uni[i in origins, j in destinations, p in products, t in periods] >= 0
    unf[i in origins, j in destinations, p in products, t in periods] >= 0
  end

  @variable(m,u[i in origins, j in destinations, p in products, t in periods, f in modes; TransportationCost[i,j,f]>0] >= 0, Int, upperbound=1e6)

  if 1 in periods
    for j in warehouses
      for p in products
        setupperbound(si[j,p,1],0)
        setlowerbound(si[j,p,1],0)
      end
    end
  end

  @constraint(m, demandsatisfaction[k in customers,p in products,t in periods],sum(x[j,k,p,t] for j in warehouses if RT[j,k,p]>0) == Dem[k,p,t])
  @constraint(m, invbalance[j in warehouses, p in products, t in periods], sf[j,p,t] == si[j,p,t] + sum(x[i,j,p,t] for i in origins if RT[i,j,p]>0) - sum(x[j,k,p,t] for k in destinations if RT[j,k,p]>0))
  @constraint(m, whupperbound[j in warehouses, p in products, t in periods], sf[j,p,t] + sum(x[j,k,p,t] for k in destinations if RT[j,k,p]>0) <= M*y[j,t,p])
  @constraint(m, plantcapacity[i in plants, t in periods, p in products], sum(x[i,j,p,t] for j in warehouses if RT[i,j,p]>0) <= capi[i,p,t])
  @constraint(m, passcapacity[i in plants, p in products, t in periods], capf[i,p,t] == capi[i,p,t] - sum(x[i,j,p,t] for j in warehouses if RT[i,j,p]>0))
  @constraint(m, transport[i in origins, j in destinations, p in products, t in periods; Route[i,j,p]>0], x[i,j,p,t] - uni[i,j,p,t]  <= sum(f*u[i,j,p,t,f] for f in modes if TransportationCost[i,j,f]>0))
  @constraint(m, leftover[i in origins, j in destinations, t in periods, p in products; Route[i,j,p]>0], unf[i,j,p,t] ==  sum(f*u[i,j,p,t,f] for f in modes if TransportationCost[i,j,f]>0) - x[i,j,p,t] - uni[i,j,p,t])
  if in(products,"PROD1")
    for i in plants
      for t in periods
        setlowerbound(capi[i,"PROD1",t],Capacity[i])
        setupperbound(capi[i,"PROD1",t],Capacity[i])
      end
    end
    for i in origins
      for j in destinations
        for t in products
          setupperbound(uni[i,j,"PROD1",t],0)
        end
      end
    end
  end

  # Linkings
  if length(products) > 1
    @constraint(m, linkp[i in plants, (k,p) in enumerate(products), t in periods; k < length(products)], capf[i,products[k],t] == capi[i,products[k+1],t])
    @constraint(m, linkp2[j in warehouses,(k,p) in enumerate(products), t in periods; k < length(products)], y[j,t,products[k]] == y[j,t,products[k+1]])
    @constraint(m, linkp3[i in origins, j in destinations, (k,p) in enumerate(products), t in periods; k < length(products)], unf[i,j,products[k],t] == uni[i,j,products[k+1],t])
  end
  if length(periods) > 1
    @constraint(m, linkt[j in warehouses, p in products, t in periods; t+1 in periods], sf[j,p,t] == si[j,p,t+1])
  end


  @expression(m, transportationCost, sum(TransportationCost[i,j,f]*u[i,j,p,t,f] for (i,j,p,t,f) in keys(u)))
  @expression(m,inventoryCost, sum(HoldingCost[j]*sf[j,p,t] for j in warehouses, p in products, t in periods))
  @expression(m,fixedCost, sum(FixedCost[j]/length(products)*y[j,t,p] for j in warehouses, p in products, t in periods))
  @expression(m, ϕ, transportationCost + inventoryCost)
  @expression(m, obj, ϕ + fixedCost)
  @objective(m, Min, obj)


  m.ext[:transportationCost] = transportationCost
  m.ext[:fixedCost] = fixedCost
  m.ext[:inventoryCost] = inventoryCost
  m.ext[:obj] = obj
  m.ext[:ϕ] = ϕ

  return m

end

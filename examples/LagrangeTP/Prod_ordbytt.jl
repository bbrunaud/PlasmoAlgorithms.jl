using Plasmo
using PlasmoAlgorithms
using Logging
using Gurobi
using DataFrames

Logging.configure(level=DEBUG)

g = PlasmoGraph()
g.solver = GurobiSolver(OutputFlag=0)

include("modelgen4.jl")

psize=20
otime=1:psize
oproducts=1:psize
###pproducts=[3, 15, 9, 2, 20, 8, 14, 12, 18, 6, 10, 4, 16, 7, 19, 13, 1, 11, 17, 5]
##firstproduct=pproducts[1]

#Create Monolitic Model
m=spmodel(oproducts,otime)
#Solve model
solve(m)
tt=getindex(m,:tt)
tt_=getvalue(tt)
#Dictionary with mean values of tt sorted by product number
tt_mean =sort( Dict((i)=> mean(tt_[s,i,t] for s in 1:length(sites) for t in 1:length(otime)) for i in 1:length(oproducts)))
#Transfor Dictionary to tuple
tt_mean_tuple = [(k,v) for (v,k) in tt_mean]
print(tt_mean_tuple)
# compute a permutation of the tuple indices that puts the tuple into sorted order by the values of tt
ord_tt=sortperm(tt_mean_tuple)
#same but with reverse order
rev_ord_tt=sortperm(tt_mean_tuple,rev=true)

#################choose between ord_tt or rev_ord_tt#######################
pord = ord_tt
#Capacity HT will be assign to firstproduct in eq3c
firstproduct=pord[1]

if 6 in psize
    function heur(mf)
        return 72827.587
    end
elseif 20 in psize
    function heur(mf)
      return 515551.12
    end
end
orders = []
gaps =[]



g = PlasmoGraph()
g.solver = GurobiSolver(OutputFlag=0)

node1 = []

for i in oproducts
    n = add_node(g)
    m = spmodel(i)
    setmodel(n,m)
    push!(node1,n)
end


node = [node1[pord[i]] for i in 1:length(pord)]
push!(orders,pord)
println("Order = $pord")

@linkconstraint(g, [s in sites, i in 1:length(oproducts)-1, t in otime],node[i][:hf][s,pord[i],t] == node[i+1][:hi][s,pord[i+1],t])


r = lagrangesolve(g, max_iterations=50,update_method=:subgradient, lagrangeheuristic=heur, initialmultipliers=:relaxation)
push!(gaps,r.gap)
#save("productorder.jld","gaps",gaps)

#end
print(firstproduct)

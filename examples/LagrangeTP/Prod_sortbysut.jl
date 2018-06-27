using Plasmo
using PlasmoAlgorithms
using Logging
using Gurobi
using DataFrames
using JLD


Logging.configure(level=DEBUG)

include("modelgen4.jl")
psize=20
otime=1:psize
oproducts=1:psize


#Dictonary with mean values of sut sorted by product number
sut_mean= sort(Dict((i) => mean(sut[s,i] for s in 1:length(sites)) for i in 1:length(oproducts)))
#Transfor Dictonary to tuple
sut_mean_tuple = [(k,v) for (v,k) in sut_mean]
print(sut_mean_tuple)
# compute a permutation of the tuple indices that puts the tuple into sorted order by the values of sut
ord_sut=sortperm(sut_mean_tuple)
#same but with reverse order
rev_ord_sut=sortperm(sut_mean_tuple,rev=true)

#choose between ord_sut or inv_ord_sut
pord = rev_ord_sut
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
save("productorder.jld","gaps",gaps)

#end
print(firstproduct)

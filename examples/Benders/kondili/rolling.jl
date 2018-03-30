using Scheduling
using Plasmo
using PlasmoAlgorithms
using Gurobi
using Ipopt
using JuMP


mipsolver = GurobiSolver(MIPGap=0.005)
nlpsolver = IpoptSolver(tol=20.0)

include("/home/bbrunaud/Scheduling/examples/kondili/kondili.jl")
n.periods = 50

g = PlasmoGraph()

t = 1

mips = []
nlps = []
node = []
bd = true

for i in 1:4
    m1 = generatemodelSTN!(n,initialperiod=t, endperiod=t+19)
    m1nl = generatemodelSTNnlp!(n, initialperiod=t, endperiod=t+19)
    @objective(m1, Min, -m1.obj)
    if bd && i < 4
        @objective(m1,Min,0)
    end
    @objective(m1nl, Min, -m1nl.obj)
    m1.solver = mipsolver
    m1nl.solver = nlpsolver
    m1.ext[:nlpmodel] = m1nl
    n1 = add_node(g,m1)
    push!(mips,m1)
    push!(nlps,m1nl)
    push!(node,n1)
    t += 10
end

@linkconstraint(g, [tk in n.tanks, s in n.materials[tk], t in 10:10:n.periods-20], node[Int(t/10)][:inv][tk,s,t] == node[Int(t/10+1)][:inv][tk,s,t])

for k in 1:length(node)-1
    add_edge(g,node[k],node[k+1])
end

#=PlasmoAlgorithms.preProcess(g)
for i in 1:3
    m = getmodel(node[i])
    θ = getindex(m, :θ)
    JuMP.fix(θ[i+1],0)
end
=#

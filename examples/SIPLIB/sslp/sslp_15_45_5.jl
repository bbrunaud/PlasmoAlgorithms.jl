using Gurobi
using JuMP

include("../generateSIPLIB.jl")

const PA = PlasmoAlgorithms
cd("examples/SIPLIB/")
g = genproblem("sslp",15,45,5)
g.solver = GurobiSolver(MIPGap=0.01)

PA.bdprepare(g)
PA.forwardstep(g,[:LP],true)

n = collect(values(getnodes(g)))
mp = n[2].model
y = getindex(mp,:y)
mp.obj

for k in keys(n[2].attributes[:cutdata])
    cd = n[2].attributes[:cutdata][k][1]
    ncd = PA.BendersCutData(cd.θk, -cd.λk, cd.xk)
    n[2].attributes[:cutdata][k][1] = ncd
end

n[2].attributes[:cutdata][1][1].λk
PA.forwardstep(g,[:LP],true)

n[2].attributes[:cutdata][6][1].λk
PA.forwardstep(g,[:LP],true)

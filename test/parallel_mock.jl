addprocs(2)
using PlasmoAlgorithms
using JuMP
using Plasmo
using Gurobi

include("fisher.jl")

n1.model.solver = GurobiSolver()
n2.model.solver = GurobiSolver()

import Base: getindex, setindex!

setindex!(n::Plasmo.PlasmoNode,value,key) = setindex!(n.attributes,value,key)
setindex!(n::Plasmo.PlasmoGraph,value,key) = setindex!(n.attributes,value,key)
getindex(n::Plasmo.PlasmoNode,key) = getindex(n.attributes,key)
getindex(n::Plasmo.PlasmoGraph,key) = getindex(n.attributes,key)


g.attributes[:levels] = Dict(1 => n1, 2 => n2)
g.attributes[:roots] = [n1]
g.attributes[:leaves] = [n2]

n2x = RemoteChannel(()->Channel{Array{Float64,1}}(1))
n1c = RemoteChannel(()->Channel{Array{PlasmoAlgorithms.CutData,1}}(1))

n2.attributes[:xin] = n2x
n1.attributes[:xout] = Dict(2 => n2x)
n1.attributes[:childvars] = Dict(2 => [xm[1], xm[2]])

put!(n2x,[100.0,99.9])

@everywhere function setx(n)
        m = n.model
        JuMP.solve(m)
        for childindex in keys(n.attributes[:childvars])
          xout = JuMP.getvalue(n.attributes[:childvars][childindex])
          old = take!(n.attributes[:xout][childindex])
          println("Taking out the trash $old")
          put!(n.attributes[:xout][childindex],xout)
          println("Putting the new in $xout")
        end
       end

n1r = remotecall_fetch(()->n1,2)
remotecall_fetch(setx,2,n1)


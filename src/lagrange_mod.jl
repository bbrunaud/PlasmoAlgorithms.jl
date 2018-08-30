
abstract type CutData end
import Base.==

struct BendersCutData <: CutData
  θk
  λk
  xk
end

struct LagrangeCutData <: CutData
  Zk
  xk
end

struct LLIntegerCutData <: CutData
  θlb
  yk
end

struct IntegerCutData <: CutData
  yk
end

(==)(cd1::BendersCutData,cd2::BendersCutData) = (cd1.θk == cd2.θk) && (cd1.λk == cd2.λk) && (cd1.xk == cd2.xk)
(==)(cd1::LLIntegerCutData,cd2::LLIntegerCutData) = (cd1.θlb == cd2.θlb) &&  (cd1.yk == cd2.yk)
(==)(cd1::IntegerCutData,cd2::IntegerCutData) = (cd1.yk == cd2.yk)


function normalizegraph(graph::ModelGraph)
    n = 1
    for node in getnodes(graph)
        m = getmodel(node)
        if m.objSense == :Max
            m.objSense = :Min
            m.obj = -m.obj
            n = -1
        end
    end
    setattribute(graph, :normalized, n)
    return n
end

graph = g
# Preprocess function
#"""
#  lgprepare(graph::PlasmoGraph)
#  Prepares the graph to apply lagrange decomposition algorithm
#"""
#function lgprepare(graph::PlasmoGraph, δ=0.5, maxnoimprove=3,cpbound=nothing)
δ=0.5
maxnoimprove=3
cpbound=1e-6
  if hasattribute(graph,:preprocessed)
    return true
  end
  n = normalizegraph(graph)
  links = getlinkconstraints(graph)
  nmult = length(links) # Number of multipliers

  setattribute(graph, :numlinks, nmult)
  setattribute(graph, :δ, δ)
  setattribute(graph, :noimprove, 0)
  setattribute(graph, :Zk, NaN)
  setattribute(graph, :maxnoimprove, maxnoimprove)
  setattribute(graph, :steptaken, false)
  setattribute(graph, :explore, [])

  # Save subproblems
  subproblems = setattribute(g, :subproblems, collect(getnodes(g)))

  # Create Lagrange Master
  graphsolver = getsolver(graph)
  ms = Model(solver=graphsolver)
  @variable(ms, η, upperbound=cpbound)
  @variable(ms, λ[1:nmult])
  @objective(ms, Max, η)

  masternode = add_node(g,ms)
  setattribute(masternode, :subλ, Dict())
  setattribute(masternode, :λout, Dict())
  setattribute(masternode, :cutdata, Dict())
  setattribute(masternode, :prevcuts, Dict())
  setattribute(graph, :masternode, masternode)

  # Each node save its initial objective and set a solver if they don't have one
  for n in subproblems
    mn = getmodel(n)
    if mn.solver == JuMP.UnsetSolver()
      mn.solver = graphsolver
    end
    mn.ext[:preobj] = mn.obj
    mn.ext[:multmap] = Dict()
    mn.ext[:varmap] = Dict()

    subindex = getindex(graph,n)
    getattribute(masternode, :subλ)[subindex] = JuMP.Variable[]

    λchannel = RemoteChannel(()->Channel{Array{Float64,1}}(1))
    getattribute(masternode, :λout)[subindex] = λchannel
    setattribute(n, :λin, λchannel)

    cutchannel = RemoteChannel(()->Channel{CutData}(10))
    getattribute(masternode, :cutdata)[subindex] = cutchannel
    getattribute(masternode, :prevcuts)[subindex] = []
    setattribute(n, :cutdataout, cutchannel)

    setattribute(n, :dualvars, JuMP.Variable[])
    setattribute(n, :dualcoeffs, Float64[])
  end

  for (i,lc) in enumerate(links)
    for j in 1:length(lc.terms.vars)
      var = lc.terms.vars[j]
      node = getnode(var)
      coeff = lc.terms.coeffs[j]
      subindex = getindex(graph, node)

      push!(getattribute(masternode,:subλ)[subindex], λ[i])
      push!(getattribute(node,:dualvars), var)
      push!(getattribute(node,:dualcoeffs), coeff)
    end
  end

  setattribute(graph, :preprocessed, true)
#end

# Multiplier Initialization
#function initialrelaxation(graph)
  if !hasattribute(graph,:mflat)
    setattribute(graph, :mflat, create_jump_graph_model(graph))
    setsolver(getattribute(graph, :mflat), getsolver(graph))
  end
  n = getattribute(graph, :normalized)
  nmult = getattribute(graph, :numlinks)
  mf = getattribute(graph, :mflat)
  JuMP.solve(mf,relaxation=true)
  masternode = getattribute(graph, :masternode)
  λ = getindex(getmodel(masternode), :λ)
  JuMP.setvalue.(λ, n*mf.linconstrDuals[end-nmult+1:end])
  return getobjectivevalue(mf)
#end

function putλ!(graph)
  masternode = getattribute(graph, :masternode)
  for idx in keys(getattribute(masternode,:subλ))
    put!(getattribute(masternode,:λout)[idx], JuMP.getvalue(getattribute(masternode,:subλ)[idx]))
  end
end
putλ!(graph)
#function solvenode(node)
node = n1
  # 1. Take λ
  m = getmodel(node)
  m.obj = m.ext[:preobj]
  λv = take!(getattribute(node,:λin))
  # Add dualized part to objective function
  coeffs = getattribute(node,:dualcoeffs)
  vars = getattribute(node,:dualvars)
  append!(m.obj.aff, AffExpr(vars,λv.*coeffs,0.0))

  # 2. Solve
  JuMP.solve(m)

  #=
  # 3. Put Zk,x
  for v in keys(m.ext[:varmap])
  val = getvalue(v)
  x[m.ext[:varmap][v]...] = val
  end

  objval = getvalue(m.ext[:lgobj])
  node.attributes[:objective] = objval
  node.attributes[:solvetime] = getsolvetime(m)

  return x, objval
end
=#

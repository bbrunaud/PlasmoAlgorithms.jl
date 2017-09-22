include("ex11Product.jl")

links = getlinkconstraints(g)
nmult = length(links)

# Assuming nodes are created in order
SP = [g.nodes[i].attributes[:model] for i in 1:length(g.nodes)]
for sp in SP
  setsolver(sp,g.solver)
end
# Capture objectives
SPObjectives = [g.nodes[i].attributes[:model].obj for i in 1:length(g.nodes)]
sense = SP[1].objSense

L = []
Obj = []
α = 2
step = 0
Zk = 0
Zkprev = 0
LB = 0

# add dualized part
function setlambda(λk)
  # Restore initial objective
  for (j,sp) in enumerate(SP)
    sp.obj = SPObjectives[j]
  end
  for l in 1:nmult
    for j in 1:length(links[l].terms.vars)
      var = links[l].terms.vars[j]
      coeff = links[l].terms.coeffs[j]
      var.m.obj += λk[l]*coeff*var
    end
  end

  map(solve,SP)
  return reduce(+,map(getobjectivevalue,SP))
end

mflat = create_flat_graph_model(g)
mflat.solver = GurobiSolver()
mflat.obj = -mflat.obj
mflat.objSense = :Max
solve(mflat,relaxation=true)

λ = mflat.linconstrDuals[end]
UB = getobjectivevalue(mflat)

objval = setlambda(λ)
push!(L,λ)
push!(Obj,objval)

λprev = λ
Zkprev = objval
dir = getvalue(links[1].terms)
λ = λprev - 2*objval/dir
objval = setlambda(λ)
push!(L,λ)
push!(Obj,objval)

if objval > Zkprev
  λ = λprev - 2/2*objval/dir
else
  λprev = λ
  Zkprev = objval
  dir = getvalue(links[1].terms)
  λ = λprev - 2*objval/dir
end

objval = setlambda(λ)
push!(L,λ)
push!(Obj,objval)

if objval > Zkprev
  λ = λprev - 2/4*objval/dir
else
  λprev = λ
  Zkprev = objval
  dir = getvalue(links[1].terms)
  λ = λprev - 2*objval/dir
end

objval = setlambda(λ)
push!(L,λ)
push!(Obj,objval)

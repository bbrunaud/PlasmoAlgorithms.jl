include("ex11Product.jl")
using Plots

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

num = 20
L1 = []
L2 = []
Obj = zeros(num,num)
ar = linspace(-20,40,num)
for (i,λ1) in enumerate(ar)
  for (j,λ2) in enumerate(ar)
    λ = (λ1,λ2)
    objval = setlambda(λ)

    Obj[i,j] = objval
  end
end

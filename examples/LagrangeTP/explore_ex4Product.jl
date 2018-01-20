using JuMP
include("ex4Product.jl")

graph = g

# Get Linkings
links = getlinkconstraints(graph)
nmult = length(links)

# Generate subproblem array
SP = [getmodel(graph.nodes[i]) for i in 1:length(graph.nodes)]
for sp in SP
  JuMP.setsolver(sp,graph.solver)
end
# Capture objectives
SPObjectives = [getmodel(graph.nodes[i]).obj for i in 1:length(graph.nodes)]
sense = SP[1].objSense

function passlambda(λk)
  # Restore initial objective
  for (j,sp) in enumerate(SP)
    sp.obj = SPObjectives[j]
  end
  # add dualized part
  for l in 1:nmult
    for j in 1:length(links[l].terms.vars)
      var = links[l].terms.vars[j]
      coeff = links[l].terms.coeffs[j]
      var.m.obj += λk[l]*coeff*var
    end
  end
end

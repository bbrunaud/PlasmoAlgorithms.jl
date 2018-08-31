
mutable struct CrossGraph
 bd
 lg
end

function CrossGraph(g::ModelGraph)
 c = CrossGraph()
 c.bd = g
 c.lg = deepcopy(g) # transformation(g)
 c
end

function crossprepare(c::CrossGraph)
  PA.bdprepare(c.bd)
  PA.lgprepare(c.lg)

  # cross-specific mappings
 end

function crosssolve(c::CrossGraph, max_iterations)
 for i in 1:max_iterations
   bendersolve(c.bd, max_iterations=1)
   benders_to_lagrange(c)
   lagrangesolve(c.lg,max_iterations=1, update_method=:cuttingplanes, callback=lagrange_to_benders)
   lagrange_to_benders(c)
 end
end

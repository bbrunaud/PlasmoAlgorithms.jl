using Plasmo 
using JuMP
using CPLEX
master = Model(solver=CplexSolver())
@variable(master, 0<=x1<=1)
@variable(master, x2, Bin)

@constraint(master, -x1 -x2 >= -1.5)

sub = Model(solver=CplexSolver()) 
@variable(sub, 0<=y1<=1)
@variable(sub, 0<=y2<=1)
@variable(sub, y3, Bin)
@variable(sub, y4, Bin)
@variable(sub, x1)
@variable(sub, x2)
# @constraint(sub, l1, x1 == 0.5)
# @constraint(sub, l2, x2==0.5)
@constraint(sub, -2*y1-3*y2-4*y3-5*y4>= -5-0.3*x1)
@constraint(sub, -6*y1-y2-3*y3-2*y4>=-10+0.3*x2)
@objective(sub, Min, -16*y1-19*y2-23*y3-28*y4)

g = PlasmoGraph()
g.solver = CplexSolver()
n1 = add_node(g)
setmodel(n1, master)
n2 = add_node(g)
setmodel(n2, sub)
@linkconstraint(g, n1[:x1]==n2[:x1])
@linkconstraint(g, n1[:x2]==n2[:x2])
edge = add_edge(g, n1, n2)
function solverootrelaxation_1norm(model)
  # sp = getmodel(node)
  sp = model
  global tmpdir = "/tmp/RootNode" 
  if length(sp.linconstrDuals) == 0
    solve(sp, relaxation=true)
  end
  lpfile = joinpath(tmpdir,"nodemodel.lp")
  writeLP(sp,lpfile)
  run(`cpxgetroot $lpfile 0 1`)
  lp = Model(solver=CPLEX.CplexSolver(CPX_PARAM_PREIND=0))
  lp.internalModel = MathProgBase.LinearQuadraticModel(lp.solver)
  # if isfile("node0.lp")
  #     MathProgBase.loadproblem!(lp.internalModel,"node0.lp")
  # else
  #   warn("Node file not found, falling back to lp relaxation")
  #   return solvelprelaxation(node)
  # end
  # # Restore Bounds
  # MathProgBase.setvarLB!(lp.internalModel,sp.colLower)
  # MathProgBase.setvarUB!(lp.internalModel,sp.colUpper)

  # MathProgBase.optimize!(lp.internalModel)

  # run(`mv node0.lp $tmpdir/`)

  # dualconstraints = node.attributes[:linkconstraints]

  # rootduals = MathProgBase.getconstrduals(lp.internalModel)
  # sp.linconstrDuals = MathProgBase.getconstrduals(lp.internalModel)[1:length(sp.linconstrDuals)]

  # λnode = getdual(dualconstraints)
  # nodebound = MathProgBase.getobjval(lp.internalModel)

  # node.attributes[:bound] = nodebound
  # node.attributes[:λ] = λnode
end

# MathProgBase.setobj!(lp.internalModel,c)

solverootrelaxation_1norm(sub)
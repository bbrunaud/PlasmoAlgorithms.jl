using Gurobi
using Plasmo
using PlasmoAlgorithms
using Logging

Logging.configure(level=DEBUG)

include("kondilimodel.jl")

# Solution of the base model
solver = GurobiSolver(OutputFlag=0)
m = fullmodel(solver)

#=solve(m)
Objective = getobjectivevalue(m)
SolveTime = getsolvetime(m)
print("BASE Model \n
Objective = $Objective \n
Time = $SolveTime")
=#
# Get First Solution
solve(m, relaxation = true)
lprelax = getobjectivevalue(m)
wr = getvalue(getindex(m,:w))
wc = Dict(k => ceil(wr[k...]) for k in keys(wr))
w1 = removeoverlaps(wc,solver)

# Get list of diverse solutions
# A is an array of w vectors for the schedule
A = [w1]
numsolutions = 20

for i in 1:numsolutions
  push!(A, diversify(A,solver))
end

# Objective value of each of the generated solutions
#fA = popfitness(A,solver)
#println("fA = $fA")

# Graph
g = PlasmoGraph()
m1 = batching(GurobiSolver(OutputFlag=0))
m2 = fullmodel(GurobiSolver(OutputFlag=0))
# Set objectives to min
@objective(m1, Min, -m1.obj)
@objective(m2, Min, -m2.obj)

n1 = add_node(g,m1)
n2 = add_node(g, m2)
add_edge(g,n1,n2)

@linkconstraint(g, [k in keys(wr)], n1[:w][k...] == n2[:w][k...])

ϵ = 10e-5
UB = Inf
LB = -Inf
max_iterations=30

#=
preProcess(g)
forwardStep(g)
for i in 1:length(A)
  initialCuts(g,[A[i]])
  backwardStep(g,:LP)
#  cutGeneration(g,n2,:Bin,θlb=lprelax)
end
=#

bendersolve(g,max_iterations=25)
#=for i in 1:max_iterations
  LB,UB = forwardStep(g)
  debug("***** ITERATION $i ***********")
  debug("*** UB = ",UB)
  debug("*** LB = ",LB)
  if abs(UB - LB)<ϵ
    debug("Converged!")
    break
  end
  backwardStep(g,:LP)
  cutGeneration(g,n2,:Bin,θlb=-2909)
end
=#
# TODO: Idea... For each solution in A, generate a cut into the master. Test with and without initialization.


## Genetic Algorithm test
#=sortpopulation(A,fA)

function newgen(A)
  NA = [zeros(length(tovec(A[1]))) for i in 1:4]
  (NA[1],NA[2]) = crossover(tovec(A[1]),tovec(A[4]))
  (NA[3],NA[4]) = crossover(tovec(A[2]),tovec(A[3]))
  for i in 1:4
    mutation(NA[i])
    todict(NA[i], A[i+4])
    A[i+4] = removeoverlaps(A[i+4],solver)
  end
end

newgen(A)
=#

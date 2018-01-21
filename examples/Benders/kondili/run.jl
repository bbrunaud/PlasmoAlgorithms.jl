using Gurobi
using Plasmo
using PlasmoAlgorithms
using Logging

Logging.configure(level=DEBUG)

include("kondilimodel.jl")

# Solution of the base model
solver = GurobiSolver()
m = fullmodel(solver)

solve(m)
Objective = getobjectivevalue(m)
SolveTime = getsolvetime(m)
print("BASE Model \n
Objective = $Objective \n
Time = $SolveTime")

# Get First Solution
solve(m, relaxation = true)
wr = getvalue(getindex(m,:w))
wc = Dict(k => ceil(wr[k...]) for k in keys(wr))
w1 = removeoverlaps(wc,solver)

# Get list of diverse solutions
# A is an array of w vectors for the schedule
A = [w1]
numsolutions = 10

for i in 1:numsolutions
  push!(A, diversify(A,solver))
end

# Objective value of each of the generated solutions
fA = map(fitness,A,[solver for i in 1:length(A)])
println("fA = $fA")

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

# TODO: Idea... For each solution in A, generate a cut into the master. Test with and without initialization.

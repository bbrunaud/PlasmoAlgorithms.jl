using JuMP

include("kondilidata.jl")

function batching(solver=JuMP.UnsetSolver())
  m = Model(solver=solver)

  @variable(m, w[i in tasks, j in UnitsForTask[i], t in periods], Bin)
  @constraint(m, onestart[j in units, t in periods], sum(w[i,j,τ] for i in TasksForUnit[j], τ in (t-ProcessingTime[i]+1):t if τ>0 && ProcessingTime[i]>0) <= 1)

  @constraint(m, d1[j in 2:3, t in periods], w[2,j,t] <= sum(w[1,1,τ] for τ in 1:t-1))
  @constraint(m, d2[j in 2:3, t in periods], w[4,j,t] <= sum(w[3,j2,τ] for τ in 1:t-2 for j2 in 2:3))
  @constraint(m, d3[j in 2:3, t in periods], w[4,j,t] <= sum(w[2,j2,τ] for τ in 1:t-2 for j2 in 2:3))
  @constraint(m, d4[t in periods], w[5,4,t] <= sum(w[4,j,τ] for τ in 1:t-1 for j in 2:3))

  @objective(m, Max, 0)
  return m
end


function fullmodel(solver)
  m = batching(solver)

  w = getindex(m,:w)
  @variable(m, b[i in tasks, j in UnitsForTask[i], t in periods] >= 0)
  @variable(m, st[s in states, t in 0:(periods[end]+1)] >= 0, upperbound=StorageCapacity[s])

  for s in states
    setlowerbound(st[s,0], InitialInventory[s])
    setupperbound(st[s,0], InitialInventory[s])
  end

    @constraint(m, unitcapacity[i in tasks, j in UnitsForTask[i], t in periods], b[i,j,t] <= UnitCapacity[j]*w[i,j,t] )

  # Produced amount
  @expression(m, produced[s in states, t in periods[1]:(periods[end]+1)], AffExpr(0.0))
  for s in states
    for t in periods[1]:(periods[end]+1)
      if length(TasksOut[s]) > 0
        for i in TasksOut[s]
          if t - ProcessingTime[i] >= periods[1]
              produced[s,t] += sum(b[i,j,t-ProcessingTime[i]]*ρout[i,s] for j in UnitsForTask[i])
          end
        end
      end
    end
  end

  # Consumed amount
  @expression(m, consumed[s in states, t in periods[1]:(periods[end]+1)], AffExpr(0.0))
  for s in states
    for t in periods
      if length(TasksIn[s]) > 0
        for i in TasksIn[s]
          consumed[s,t] += sum(b[i,j,t]*ρin[i,s] for j in UnitsForTask[i])
        end
      end
    end
  end

  @constraint(m, massbalance[s in states, t in periods[1]:(periods[end]+1)], st[s,t] == st[s,t-1] + produced[s,t] - consumed[s,t] )

  @objective(m, Max, sum(Price[s]*st[s,periods[end]+1] for s in states))

  return m
end

function diversemodel(A::Array{<:Any,1},solver)
    m = batching(solver)
    w = getindex(m,:w)

    @objective(m, Max, sum(w[i...] for wk in A for i in keys(w) if wk[i...] == 0)
    + sum(1-w[i...] for wk in A for i in keys(w) if wk[i...] == 1))

    return m
end

function diversify(A, solver)
    m = diversemodel(A, solver)
    solve(m)

    return getvalue(getindex(m,:w))
end


function removeoverlaps(wr,solver)
    m = diversemodel([wr],solver)
    m.obj = -m.obj

    solve(m)

    return getvalue(getindex(m,:w))
end

function fitness(wv,solver)
    m = fullmodel(solver)
    w = getindex(m,:w)
    for k in keys(w)
        JuMP.fix(w[k...],wv[k...])
    end

    status = solve(m)
    if status == :Optimal
        return getobjectivevalue(m)
    else
        return -1e5
    end
end

"""
bendersolve
"""
function bendersolve(graph::Plasmo.PlasmoGraph; max_iterations::Int64=10, cuts::Array{Symbol,1}=[:LP], ϵ=1e-5,UB=Inf,UBupdatefrequency=1)
  starttime = time()
  s = Solution(method=:benders)

  bdprepare(graph,UBupdatefrequency)
  n = graph.attributes[:normalized]

  solve(graph.attributes[:mflat])
  LB = n*getobjectivevalue(graph.attributes[:mflat])
  UB = n*UB

  for i in 1:max_iterations
    tic()
    LB,UB = forwardstep(graph::PlasmoGraph, cuts::Array{Symbol,1})

    if abs(UB-LB) < ϵ
      s.termination = "Optimal"
      return s
    end

    # Check time limit
    tstamp = starttime - time()
    if tstamp > timelimit
      s.termination = "Time Limit"
      return s
    end
    saveiteration(s,tstamp,[n*LB,n*UB,toc(),tstamp],n)
    printiterationsummary(s,singleline=false)
  end

  s.termination = "Max Iterations"
  return s
end

function forwardstep(graph::PlasmoGraph, cuts::Array{Symbol,1})
  levels = graph.attributes[:levels]
  numlevels = length(levels)
  for level in 1:numlevels
    nodeslevel = levels[level]
    for node in nodeslevel
      solvenode(node,cuts)
    end
  end
end

function solvenode(node::PlasmoNode, cuts::Array{Symbol,1})
  # 1. Add cuts
  # 2. Take x
  # 3. solve
  # 4. put x
  # 5. put λ
end

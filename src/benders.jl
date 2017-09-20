function benderssolve(graph::PlasmoGraph;max_iterations=2)
# bendersolve(g, max_iterations=10)

  ########## 1. Initialize ########
  tic()
  starttime = time()
  # Results outputs
  df = DataFrame(Iter=[],Time=[],α=[],step=[],UB=[],LB=[],Hk=[],Zk=[],Gap=[])
  res = Dict()

  
  push!(df,[iter,round(time()-starttime),α,step,UB,LB,Hk,Zk,gap])
  return res, df
end

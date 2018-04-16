###TODO ask braulio what way we should define models
###Might be easier to define τ, cT, dT, H and then auto build graph
###or take in original model and break it down (by indetifying complicating variables)
###right now cuts depend on knowing τ, cT, dT, H explicitly



using JuMP
using Plasmo

function fixBSP(graph, BSPs, xbar)
    childLinks = graph.attributes[:childlinks]
    for s = 1:length(BSPs)
        BSP = BSPs[s]
        BSPmodel = getmodel(BSP)
        xs = []
        for childLink in childLinks[BSP]
          x = getindex(BSPmodel,:valbar)
          push!(xs,x[childLink])
        end
        for i = 1:length(xs)
            setupperbound(xs[i],xbar[i])
            setlowerbound(xs[i],xbar[i])
        end
    end
end

function fixLSP(LSPs, origObjs, H, λ)
    for s = 1:length(LSPs)
        LSP = LSPs[s]
        origObj = origObjs[s]
        x = getindex(LSP, :x)
        @objective(LSP, Min, origObj + sum(x.*H[s]*λ[s]))
    end
    return λ
end

function isChildNode(g::Plasmo.PlasmoGraph, n1::PlasmoNode, n2::PlasmoNode)
  ##Checks if n1 is a child node of n2
  for node in LightGraphs.out_neighbors(g.graph,getindex(g,n2))
    if (n1 == g.nodes[node]) return true
      return true
    end
  end
  return false
end

function numChildNodes(g::PlasmoGraph, n1::PlasmoNode)
  return length(LightGraphs.out_neighbors(g.graph,getindex(g,n1)))
end

function numParentNodes(g::PlasmoGraph, n1::PlasmoNode)
  return length(LightGraphs.in_neighbors(g.graph,getindex(g,n1)))
end

function bendersPrep(graph::PlasmoGraph)
  links = getlinkconstraints(graph)
  numLinks = length(links)
  numNodes = length(graph.nodes)

  graph.attributes[:links] = Dict()
  graph.attributes[:childlinks] = Dict()
  graph.attributes[:duals] = Dict()
  dualMap = graph.attributes[:duals]
  linkIndex = graph.attributes[:links]
  childLinkIndex = graph.attributes[:childlinks]

  #Add each node as a key to the dictionary
  for index in 1:length(graph.nodes)
    node = graph.nodes[index]
    nodeModel = getmodel(node)
    nodeModel.ext[:preobj] = nodeModel.obj
    linkIndex[node] = []
    childLinkIndex[node] = []
    #Add valbar to child nodes
    if numParentNodes(graph,node) != 0
      sp = getmodel(node)
      @constraintref dual[1:numLinks]
      @variable(sp, valbar[1:numLinks])
    end
    #Add theta to parent nodes
    if numChildNodes(graph,node) != 0
      mp = getmodel(node)
      @variable(mp,θ[1:numNodes] >= -1e6)
      for node in LightGraphs.out_neighbors(graph.graph,getindex(graph,node))
        childNode = graph.nodes[node]
        childNodeIndex = getindex(graph,childNode)
        θ = getindex(mp,:θ)
        mp.obj += θ[childNodeIndex]
      end
    end
  end
  for link in 1:numLinks
    #Take the two variables of the constraint
    var1 = links[link].terms.vars[1]
    var2 = links[link].terms.vars[2]
    #Determine which nodes they belong to
    nodeV1 = getnode(var1)
    nodeV2 = getnode(var2)
    #Set the order of the nodes
    if isChildNode(graph,nodeV1,nodeV2)
      childNode = nodeV1
      childvar = var1
      parentNode = nodeV2
    elseif isChildNode(graph,nodeV2,nodeV1)
      childNode = nodeV2
      childvar = var2
      parentNode = nodeV1
    end
    #Add theta[childNode] to the parent node objective
    mp = getmodel(parentNode)
    childNodeIndex = getindex(graph,childNode)
    # θ = getindex(mp,:θ)
    # mp.obj += θ[childNodeIndex]
    #Add linking constraint to the parent node dictionary
    linkList = linkIndex[parentNode]
    childLinkList = childLinkIndex[childNode]
    push!(linkList,link)
    push!(childLinkList,link)
  end
  for child in keys(childLinkIndex)
    sp = getmodel(child)
    childLinkList = childLinkIndex[child]
    @constraintref dual[1:length(childLinkList)]
    for i in 1:length(childLinkList)
      link = childLinkList[i]
      var1 = links[link].terms.vars[1]
      var2 = links[link].terms.vars[2]
      #Determine which nodes they belong to
      nodeV1 = getnode(var1)
      nodeV2 = getnode(var2)
      #Set the order of the nodes
      if isChildNode(graph,nodeV1,nodeV2)
        childNode = nodeV1
        childvar = var1
        parentNode = nodeV2
      elseif isChildNode(graph,nodeV2,nodeV1)
        childNode = nodeV2
        childvar = var2
        parentNode = nodeV1
      end
      valbar = getindex(sp,:valbar)
      dual[i] = @constraint(sp,valbar[link] - childvar == 0)
    end
    dualMap[child] = dual
  end
  graph.attributes[:BMP] = []
  graph.attributes[:BSPs] = []
  #Create dictionary to keep track of levels of nodes
  graph.attributes[:levels] = Dict()
  #Iterate through every node to check for root/leaf nodes
  for nodeIndex in 1:length(graph.nodes)
    node = graph.nodes[nodeIndex]
    #If the node does not have parents it is the BMP
    if numParentNodes(graph,node)==0
      push!(graph.attributes[:BMP],node)
      #Root nodes are the first level
    end
    #If the node does not have children it is a BSPs
    if numChildNodes(graph,node)==0
      push!(graph.attributes[:BSPs],node)
    end
  end
end

function lagrangePrep(graph::PlasmoGraph, τ)
  links = getlinkconstraints(graph)
  numLinks = length(links)
  numNodes = length(graph.nodes)

  graph.attributes[:links] = Dict()
  graph.attributes[:childlinks] = Dict()
  graph.attributes[:duals] = Dict()
  dualMap = graph.attributes[:duals]
  linkIndex = graph.attributes[:links]
  childLinkIndex = graph.attributes[:childlinks]

  #Add each node as a key to the dictionary
  for index in 1:length(graph.nodes)
    node = graph.nodes[index]
    nodeModel = getmodel(node)
    nodeModel.ext[:preobj] = nodeModel.obj
    linkIndex[node] = []
    childLinkIndex[node] = []
    #Add valbar to child nodes
    if numParentNodes(graph,node) != 0
      sp = getmodel(node)
      @constraintref dual[1:numLinks]
      @variable(sp, valbar[1:numLinks])
    end
  end
  for link in 1:numLinks
    #Take the two variables of the constraint
    var1 = links[link].terms.vars[1]
    var2 = links[link].terms.vars[2]
    #Determine which nodes they belong to
    nodeV1 = getnode(var1)
    nodeV2 = getnode(var2)
    #Set the order of the nodes
    if isChildNode(graph,nodeV1,nodeV2)
      childNode = nodeV1
      childvar = var1
      parentNode = nodeV2
    elseif isChildNode(graph,nodeV2,nodeV1)
      childNode = nodeV2
      childvar = var2
      parentNode = nodeV1
    end
    #Add theta[childNode] to the parent node objective
    mp = getmodel(parentNode)
    childNodeIndex = getindex(graph,childNode)
    # θ = getindex(mp,:θ)
    # mp.obj += θ[childNodeIndex]
    #Add linking constraint to the parent node dictionary
    linkList = linkIndex[parentNode]
    childLinkList = childLinkIndex[childNode]
    push!(linkList,link)
    push!(childLinkList,link)
  end
  for child in keys(childLinkIndex)
    sp = getmodel(child)
    childLinkList = childLinkIndex[child]
    @constraintref dual[1:length(childLinkList)]
    for i in 1:length(childLinkList)
      link = childLinkList[i]
      var1 = links[link].terms.vars[1]
      var2 = links[link].terms.vars[2]
      #Determine which nodes they belong to
      nodeV1 = getnode(var1)
      nodeV2 = getnode(var2)
      #Set the order of the nodes
      if isChildNode(graph,nodeV1,nodeV2)
        childNode = nodeV1
        childvar = var1
        parentNode = nodeV2
      elseif isChildNode(graph,nodeV2,nodeV1)
        childNode = nodeV2
        childvar = var2
        parentNode = nodeV1
      end
      valbar = getindex(sp,:valbar)
      dual[i] = @constraint(sp,valbar[link] - childvar == 0)
    end
    dualMap[child] = dual
  end
  graph.attributes[:LMP] = []
  graph.attributes[:LSPs] = []
  #Iterate through every node to check for root/leaf nodes
  for nodeIndex in 1:length(graph.nodes)
    node = graph.nodes[nodeIndex]
    #If the node does not have parents it is the BMP
    if numParentNodes(graph,node)==0
      push!(graph.attributes[:LMP],node)
      #Root nodes are the first level
    end
    #If the node does not have children it is a BSPs
    if numChildNodes(graph,node)==0
      push!(graph.attributes[:LSPs],node)
    end
  end
  numLSPs = length(graph.attributes[:LSPs])
  origModel = getmodel(graph.attributes[:LMP][1])
  origObj = origModel.obj
  numNodes = length(graph.nodes)
  LMPmodel = Model(solver = GurobiSolver())
  @variable(LMPmodel, κ[1:numLSPs]>=0)
  @variable(LMPmodel, η)
  @constraint(LMPmodel, η<=sum(κ))
  @objective(LMPmodel, Max, η)
  setmodel(graph.nodes[1], LMPmodel)
  println("****SOS*****")
  for nodeIndex = 2:numNodes
        m = getmodel(graph.nodes[nodeIndex])
        println("********", nodeIndex)
        println(1/(numNodes-1)*origObj)
        m.obj += 1/(numNodes-1)*origObj
        print(m)
  end
  return graph
end

function preProcess(graph::PlasmoGraph)
  lGraph = deepcopy(graph)
  if get(graph.attributes,:preprocessed,false)
    return 0
  end
  bendersPrep(graph)
  lGraph = lagrangePrep(lGraph)
  BMP, BSPs = graph.attributes[:BMP][1], graph.attributes[:BSPs]
  return graph, lGraph
end

function initialize(BSPs, LSPs, H)
    origObjs = []
    fixBSP(graph, BSPs, [0,0])
    for s = 1:length(LSPs)
        push!(origObjs, LSPs[s].obj)
    end
    fixLSP(LSPs, origObjs, H, [0,0])
    return origObjs, [0,0]
end

function BSPsolve(graph, BSPs)
    dualMap = graph.attributes[:duals]

    xhat = []
    zk = []
    μ = []
    for s = 1:length(BSPs)
        BSP = BSPs[s]
        BSPmodel = getmodel(BSP)
        solve(BSPmodel)

        valbar = getvalue(getindex(BSPmodel, :valbar))
        push!(xhat, valbar)

        z = getobjectivevalue(BSPmodel)
        push!(zk, z)

        dualCon = dualMap[BSP]
        μs = getdual(dualCon)

        push!(μ, μs)

    end
    return μ, zk, xhat
end

function LSPsolve(LSPs)
    xtilde = []
    ytilde = []
    for s = 1:length(BSPs)
        LSP = LSPs[s]
        solve(LSP)

        x = getindex(LSP, :x)
        push!(xtilde, getvalue(x))

        y = getindex(LSP, :y)
        push!(ytilde, getvalue(y))
    end
    return xtilde, ytilde
end

function BSPtoLMP(LMP::JuMP.Model, zk, H, xhat)
    κ = getindex(LMP, :κ)
    λ = getindex(LMP, :λ)
    for s = 1:length(zk)
        @constraint(LMP, κ[s] <= zk[s]+sum(λ[s].*H[s].*xhat[s]))
    end
end

function BSPtoBMP(graph, BSPs, BMP, μ)
    childLinks = graph.attributes[:childlinks]
    BMPmodel = getmodel(BMP)
    θ = getindex(BMPmodel, :θ)
    x = getindex(BMPmodel, :x)
    for s = 1:length(BSPs)
        BSP = BSPs[s]
        BSPmodel = getmodel(BSP)
        valbars = []
        for childLink in childLinks[BSP]
          valbar = getindex(BSPmodel,:valbar)
          push!(valbars,getvalue(valbar[childLink]))
        end
        @constraint(BMPmodel, θ[s] >= getobjectivevalue(BSPmodel)- sum(μ[s].*valbars)+sum(μ[s].*x))

    end
end

function LSPtoLMP(LMP, τ, cT, dT, H, xtilde, ytilde)
    κ = getindex(LMP, :κ)
    λ = getindex(LMP, :λ)

    for s = 1:length(xtilde)
        @constraint(LMP, κ[s] <= sum(τ[s]*cT.*xtilde[s])+sum(τ[s]*dT[s].*ytilde[s])+λ[s]*sum(H[s].*xtilde[s]))
    end
end

function LSPtoBMP(LSPs, BMP, H, λ)
    θ = getindex(BMP, :θ)
    x = getindex(BMP, :x)
    for s = 1:length(LSPs)
        LSP = LSPs[s]
        @constraint(BMP, θ[s] >= getobjectivevalue(LSP)-sum(λ[s].*H[s].*x))
    end

end

function LMPtoLSP(LMP, LSPs, origObjs, H)
    solve(LMP)
    λ = getvalue(getindex(LMP, :λ))
    fixLSP(LSPs, origObjs, H, λ)
    return λ
end



function BMPtoBSP(graph, BMP, BSPs)
    BMPmodel  = getmodel(BMP)
    print(BMPmodel)
    solve(BMPmodel)
    x = getindex(BMPmodel, :x)
    fixBSP(graph, BSPs, getvalue(x))
end

function crossSolve(graph::PlasmoGraph, max_iterations::Int64)
  return preProcess(graph)
    #BMP, BSPs = preProcess(graph)
    for i = 1:max_iterations
        #μ, zk, xhat = BSPsolve(graph,BSPs)
        #BSPtoBMP(graph, BSPs, BMP, μ)
        #BSPtoLMP(LMP, zk, H, xhat)
        #xtilde, ytilde = LSPsolve(LSPs)
        #LSPtoLMP(LMP, τ, cT, dT, H, xtilde, ytilde)
        #LSPtoBMP(LSPs, BMP, H, λ)
        #BMPtoBSP(graph, BMP, BSPs)
        #print("h8")
        #λ = LMPtoLSP(LMP, LSPs, origObjs, H)
        #print("h9")
    end
end

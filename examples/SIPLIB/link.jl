using Plasmo

# Collect the nodes from the original graph
n = collect(values(getnodes(g)))

g2 = PlasmoGraph()
n1 = add_node(g2,getnodes(g)[2].model)
n2 = add_node(g2,getnodes(g)[1].model)
y = getindex(n1.model,:y)

@linkconstraint(g2, [i in keys(y)], n1[:y][i] == n2[:x][i])




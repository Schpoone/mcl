from mcl_clustering import networkx_mcl
import networkx as nx
import matplotlib.pyplot as plt

G = nx.Graph()
G.add_node(1, num="one")
G.add_node(2, num="two")
G.add_node(3, num="three")
G.add_node(4, num="four")
G.add_node(5, num="five")
G.add_node(6, num="six")
G.add_node(7, num="seven")
G.add_node(8, num="eight")
G.add_node(9, num="nine")
G.add_node(0, num="zero")
G.add_edge(1,2)
G.add_edge(1,3)
G.add_edge(1,4)
G.add_edge(3,4)
G.add_edge(2,3)
G.add_edge(2,0)
G.add_edge(2,5)
G.add_edge(0,8)
G.add_edge(0,6)
G.add_edge(5,8)
G.add_edge(4,5)
G.add_edge(5,6)
G.add_edge(8,6)
G.add_edge(5,9)
G.add_edge(6,7)
G.add_edge(5,7)
G.add_edge(6,9)
G.add_edge(7,9)

M, clusters = networkx_mcl(G, inflate_factor = 1.8,
                   max_loop = 100)

print clusters
print nx.get_node_attributes(G, "num")

#nx.draw_networkx(G)
#plt.show()

#    M = output matrix
#    clusters = dict with keys = [<cluster id>] values = [<vertex id>]

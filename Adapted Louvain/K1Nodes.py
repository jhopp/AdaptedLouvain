from uuid import getnode
import community as community_louvain
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities, modularity
import os
from timeit import default_timer as timer

# Returns nodes of degree 1
def getNodesK1(graph):
    k1Nodes = 0
    for node, degree in dict(graph.degree()).items():
        if degree == 1:
            k1Nodes+=1
    return k1Nodes

# Returns network generated from .txt file
def generateNetworkFromFile(network_name):
    directory = os.getcwd()
    path_to_networks = directory + "\\..\\Networks\\"
    if(network_name.startswith("fb")):
        G = nx.read_edgelist(path_to_networks + network_name, delimiter=' ', create_using=nx.Graph())
    else:
        G = nx.read_edgelist(path_to_networks + network_name, delimiter='\t', create_using=nx.Graph())
    print("initialized network: " + network_name)
    return G

network_name = "fb107" + ".edgelist.txt"
G = generateNetworkFromFile(network_name)
n=G.number_of_nodes()
print(2*G.number_of_edges()/G.number_of_nodes())
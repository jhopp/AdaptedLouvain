import community as community_louvain
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import networkx as nx
import os
from timeit import default_timer as timer
import random as rnd
import pandas as pd
from scipy.stats import powerlaw, ttest_ind
def main():
    # Parameters
    runs = 20
    avgK = 15
    fractions = [0.05, 0.2, 0.5, 0.8]
    nodes = [100, 1000, 5000, 10000]
    debug = False
    #

    col = ['Fraction'] + nodes
    table = pd.DataFrame(columns=col)

    for i in range(0,len(fractions)):
        standardTime = {}
        adaptedTime = {}
        standardMod = {}
        adaptedMod = {}
        standardTime['Fraction'] = fractions[i]
        adaptedTime['Fraction'] = fractions[i]
        standardMod['Fraction'] = fractions[i]
        adaptedMod['Fraction'] = fractions[i]
        for j in range(len(nodes)):
            degreeDist = getDegreeDistribution(nodes[j], fractions[i], avgK, rand_seed=j)
            g = generateConfigurationNetwork(degreeDist, rand_seed = j)
            #g = generateRandomNetwork()

            # Print network info
            print("maximum degree :" + str(max(degreeDist)))
            print("number of nodes: {}".format(g.number_of_nodes()))
            print("number of edges: {}".format(g.number_of_edges()))
            print("average degree:" + str(2*g.number_of_edges()/g.number_of_nodes()))
            
            # Accumulator variables
            averageTimeNormalAcc = 0
            averageTimeAdaptedAcc = 0
            averageModNormalAcc = 0
            averageModAdaptedAcc = 0

            tModularities = []
            tAdaptModularities =[]
            for ii in range(0, runs):
                # Louvain
                time1, partition1 = runAndTime(community_louvain.best_partition, g, random_state=ii)
                modularity1 = community_louvain.modularity(partition1, g)

                # Adapted Louvain
                x = g.copy()
                time2, (partition2, k1Nodes) = runAndTime(adaptedLouvain, x, rand_state=ii)
                modularity2 = community_louvain.modularity(partition2, g)

                # Print stats
                if (debug):
                    print("Number of k1 nodes: " + str(len(k1Nodes)))
                    print("Louvain modularity:" + str(modularity1))
                    print("Adapted Louvain modularity: " + str(modularity2))
                    print("Louvain runtime: " + str(time1))  
                    print("Adapted Louvain runtime: " + str(time2))
        
                if(ii%5 ==0):
                    print(ii)

                # Update accumulators
                averageTimeNormalAcc += time1
                averageTimeAdaptedAcc += time2
                averageModNormalAcc += modularity1
                averageModAdaptedAcc += modularity2
                #
                tModularities.append(modularity1)
                tAdaptModularities.append(modularity2)

            # Print results
            print("Number of k1 nodes: " + str(len(k1Nodes)))
            print("Average time standard louvain: " + str(averageTimeNormalAcc/runs))
            print("Average time adapted louvain: " + str(averageTimeAdaptedAcc/runs))
            print("Average mod standard louvain: " + str(averageModNormalAcc/runs))
            print("Average mod adapted louvain: " + str(averageModAdaptedAcc/runs))
            standardTime[nodes[j]] = round(averageTimeNormalAcc/runs, 4)
            adaptedTime[nodes[j]] = round(averageTimeAdaptedAcc/runs, 4)
            standardMod[nodes[j]] = round(averageModNormalAcc/runs, 4)
            adaptedMod[nodes[j]] = round(averageModAdaptedAcc/runs, 4)
            tstat, pval = ttest_ind(tModularities, tAdaptModularities)
            print("statistic: " + str(tstat) + " p-value: " + str(pval))
            #

        # Add results to table
        table = table.append(standardTime, ignore_index = True)
        table = table.append(adaptedTime, ignore_index = True)
        table = table.append(standardMod, ignore_index = True)
        table = table.append(adaptedMod, ignore_index = True)
        #
    table.to_csv("output.csv", sep = ";", decimal=",", index = False)

# Run function and return runtime and result
def runAndTime(func, *args, **kwargs):
    startTime = timer()
    res = func(*args, **kwargs)
    stopTime = timer()
    return stopTime - startTime, res

# Get nodes of degree 1 from graph with its neighbours
def getNodesK1(graph):
    k1Nodes = []
    neighbours = []
    for node, degree in dict(graph.degree()).items():
        if degree == 1:
            k1Nodes.append(node)
            neighbours.append(list(graph[node])[0])
    return k1Nodes, neighbours

# Run adapted Louvain algorithm (without k1) and return resulting partition
def adaptedLouvain(graph, rand_state=None):
    k1Nodes, neighbours = getNodesK1(graph)
    graph.remove_nodes_from(k1Nodes)
    partition = community_louvain.best_partition(graph,random_state= rand_state)
    partition = addK1ToPartition(k1Nodes,neighbours, partition)
    return partition, k1Nodes

# Add nodes of degree 1 to partition, creating new communities and return new partition
def addK1ToPartition(k1Nodes,neighbours, partition):
    community_counter = max(partition.values())
    for i in range(0, len(k1Nodes)):
        key = neighbours[i]
        if key in partition:
            partition[k1Nodes[i]] = partition[neighbours[i]]
        else:
            community_counter += 1
            partition[k1Nodes[i]] = community_counter
    return partition

# Returns network generated from file
def generateNetworkFromFile(network_name):
    directory = os.getcwd()
    path_to_networks = directory + "\\..\\Networks\\"
    if(network_name.startswith("fb")):
        g = nx.read_edgelist(path_to_networks + network_name, delimiter=' ', create_using=nx.Graph())
    else:
        g = nx.read_edgelist(path_to_networks + network_name, delimiter='\t', create_using=nx.Graph())
    print("initialized network: " + network_name)
    return g

# Returns degree distribution with a fraction of nodes being degree 1
def getDegreeDistribution(nodes, fraction, averageDegree, rand_seed = None):
    totalEdges = nodes * averageDegree / 2
    k1Nodes = fraction * nodes
    remainingNodes = nodes * (1 - fraction)
    newAverageDegree = ((totalEdges - 0.5*k1Nodes) / remainingNodes) * 2 - 1
    newGamma = (1 - 2*newAverageDegree) / (1-newAverageDegree)
    
    graphAvgDegree = 0
    attempts = 0
    bestAvgDegreeDiff = 99999
    z = nx.utils.powerlaw_sequence(int(remainingNodes), newGamma, seed =rand_seed)
    averageAvgDegree = 0

    # attempt to find good distribution
    while(bestAvgDegreeDiff > 0.2 and attempts < 10000):
        newGraph = nx.utils.powerlaw_sequence(int(remainingNodes), newGamma, seed = rand_seed)
        graphAvgDegree = sum(newGraph)/len(newGraph)
        avgDegreeDiff = abs(newAverageDegree - graphAvgDegree)

        if avgDegreeDiff < bestAvgDegreeDiff:
            bestAvgDegreeDiff = avgDegreeDiff
            z = newGraph.copy()
        attempts += 1
        averageAvgDegree += avgDegreeDiff
    print(averageAvgDegree/attempts)
    if(attempts >= 10000):
        print("failed to find good distribution, using best found")
    print("attempts:" + str(attempts))

    # increase degrees by 1
    for i in range(0,len(z)):
        z[i]= z[i] + 1

    # round numbers
    for i in range(0,len(z)):
        r = rnd.random()
        if r > z[i] - int(z[i]):
            z[i] = int(z[i])
        else:
            z[i]= int(z[i] + 1)

    # add nodes degree 1
    for i in range(0,int(k1Nodes)):
        z.append(1)
    if(sum(z)%2 != 0):
        z[0]+=1
    return z

def generateConfigurationNetwork(degreeDist, rand_seed = None):
    return nx.configuration_model(degreeDist, seed= rand_seed)

if __name__ == '__main__':
    main()

# Returns randomly generated network
def generateRandomNetwork():
    G = nx.fast_gnp_random_graph(10000,1/20000.0)
    print("initialized random network")
    return G
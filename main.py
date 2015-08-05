from mcl_clustering import networkx_mcl
import networkx as nx
import matplotlib.pyplot as plt
import random
import time

"""
Reads data from a file to create a graph structure

@param fileName, the name of the file to be read
@return a graph created from the data in the file
"""
def readDataToGraph(fileName):
    file = open(fileName)
    #2-d array where the columns are interaction type, gene 1, and gene 2
    rows = []
    for line in file:
        cols = line.split("\t")
        rows.append(cols)
    file.close()
    graph = nx.Graph()
    for row in rows:
        interaction = row[0]
        gene1 = row[1]
        gene2 = row[2]
        if interaction == '""' or gene1 == '""' or gene2 == '""':
            continue
        #The block of "if not" statements makes sure there are
        #no duplicate nodes or edges
        if not graph.has_node(gene1):
            graph.add_node(gene1)
        if not graph.has_node(gene2):
            graph.add_node(gene2)
        if not graph.has_edge(gene1, gene2):
            graph.add_edge(gene1, gene2)
    return graph

"""
Takes a simple random sample of interaction data to create a graph

@param fileName, the name of the file to be read
@param sampleSize, the number of rows in the sample
@return a graph created from the sample of data in the file
"""
def sampleData(fileName, sampleSize):
    file = open(fileName)
    #2-d array where the columns are interaction type, gene 1, and gene 2
    rows = []
    for line in file:
        cols = line.split("\t")
        rows.append(cols)
    file.close()
    graph = nx.Graph()
    counter = sampleSize
    nodeID = 0
    while counter > 0:
        previousNodes = graph.nodes()
        previousEdges = graph.edges()
        randomRowNumber = int(random.random()*len(rows))
        row = rows[randomRowNumber]
        interaction = row[0].strip()
        gene1 = row[1].strip()
        gene2 = row[2].strip()
        if interaction == '""' or gene1 == '""' or gene2 == '""':
            continue
        #The block of "if not" statements makes sure there are
        #no duplicate nodes or edges
        genesByNodeID = nx.get_node_attributes(graph, "gene") #{nodeID:geneIdentifier,...}
        if gene1 not in genesByNodeID.values():
            graph.add_node(nodeID, gene=gene1)
            nodeID1 = nodeID
            nodeID+=1
        if gene2 not in genesByNodeID.values():
            graph.add_node(nodeID, gene=gene2)
            nodeID2 = nodeID
            nodeID+=1
        if not graph.has_edge(nodeID1, nodeID2):
            graph.add_edge(nodeID1, nodeID2)
        if graph.nodes() == previousNodes and graph.edges() == previousEdges:
            continue
        counter-=1
    return graph

"""
Generates a graph with random edges between nodes given a pre-generated graph. The nodes
of the original graph will be kept but any edges will be removed and replaced with
the same number of random ones.

@param graph, the original graph
@return a graph with random edges
"""
def generateRandomGraph(graph):
    numEdges = graph.size()
    nodes = graph.nodes()
    edges = graph.edges()
    for edge in edges:
        graph.remove_edge(*edge)
    currentEdges = 0
    while currentEdges < numEdges:
        node1 = nodes[int(random.random()*len(nodes))]
        node2 = nodes[int(random.random()*len(nodes))]
        if node1 == node2:
            continue
        else:
            graph.add_edge(node1, node2)
            currentEdges+=1
    return graph

"""
Generates a complete graph given a pre-generated graph. A complete graph is one where each
node has an edge between every other node. The nodes of the original graph will be kept.
"""
def generateCompleteGraph(graph):
    nodes = graph.nodes()
    for node1 in nodes:
        for node2 in nodes:
            if not graph.has_edge(node1, node2) and node1 != node2:
                graph.add_edge(node1, node2)
    return graph

"""
Display some basic statistics about a graph

@param graph
"""
def displayGraphStats(graph):
    print "Number of nodes:", graph.order()
    print "Number of edges:", graph.size()
    print "Number of nodes with self loops:", len(graph.nodes_with_selfloops())
    print "Average Degree:", sum(graph.degree().values())/float(graph.order())

"""
Converts a dictionary of clusters with node IDs to a dictionary of clusters
with gene identifiers.

Necessary to convert the output of the MCL algorithm to an in-context output

@param graph, the original graph which relates node IDs to gene identifiers
@param clusters, the dictionary of clusters to be converted
@return a dictionary of clusters using gene identifiers
"""
def convertClustersOfNodeIDsToClustersOfGeneIdentifiers(graph, clusters): #{clusterIdentifier:nodeIDs,...}
    clustersWithGeneIdentifiers = {} #{clusterIdentifier:geneIdentifiers,...}
    genesByNodeID = nx.get_node_attributes(graph, "gene") #{nodeID:geneIdentifier,...}
    for clusterIdentifier in clusters:
        cluster = clusters[clusterIdentifier]
        clusterWithGeneIdentifiers = []
        for nodeID in cluster:
            geneIdentifier = genesByNodeID[nodeID]
            clusterWithGeneIdentifiers.append(geneIdentifier)
        clustersWithGeneIdentifiers[clusterIdentifier] = clusterWithGeneIdentifiers
    return clustersWithGeneIdentifiers

"""
Removes duplicate clusters in a dictionary of clusters

@param clusters, the dictionary of clusters with duplicate clusters
@return the dictionary of clusters without duplicate clusters
"""
def removeDuplicateClusters(clusters):
    noDupes = {}
    clusterCounter = 0
    for cluster in clusters.values():
        for gene in cluster:
            if cluster in noDupes.values():
                continue
            noDupes[clusterCounter] = cluster
            clusterCounter+=1
    return noDupes

"""
Display some basic statistics about a set of clusters

@param clusters, the list of clusters
"""
def displayClusterStats(clusters):
    print "Number of clusters:", len(clusters)
    i = 0
    sumClusterSize = 0
    for cluster in clusters.values():
        sumClusterSize+=len(cluster)
        if len(cluster) == 1:
            i+=1
    print "Average cluster size:", sumClusterSize/float(len(clusters))
    print "Number of lone clusters:", i

"""
Read pathway data so that it can be compared with
the output clusters.

@return the data in a dictionary
"""
def readPathwayData(fileName):
    file = open(fileName)
    #2-d array where the columns are pathway identifier, gene
    rows = []
    for line in file:
        cols = line.split("\t")
        rows.append(cols)
    file.close()
    del rows[0] #Those are headers
    pathways = {}
    for row in rows:
        pathID = row[0]
        gene = row[1]
        if pathID in pathways:
            pathways[pathID].append(gene)
        else:
            pathways[pathID] = [gene]
    return pathways

"""
Display some basic statistics about a set of pathways.

@param pathways, the dict of pathways
"""
def displayPathwayStats(pathways):
    print "Number of Pathways:", len(pathways)
    sumPathwaySize = 0
    for pathway in pathways.values():
        sumPathwaySize+=len(pathway)
    print "Total Number of Genes:", sumPathwaySize
    print "Average Pathway Size:", sumPathwaySize/float(len(pathways))

"""
Helper function to make sure cluster information can be obtained from each individual gene
Will return two lists and each index has the gene in one list and the cluster identifier in the other list

@param clusters, the dictionary of clusters
@return an ordered list of genes and an ordered list of cluster identifiers
"""
def assignClusterToGenes(clusters):
    allGenes = []
    clusterIdentifiers = []
    for clusterIdentifier in clusters:
        for gene in clusters[clusterIdentifier]:
            allGenes.append(gene)
            clusterIdentifiers.append(clusterIdentifier)
    return allGenes, clusterIdentifiers

"""
Helper function to make sure pathway information can be obtained from each individual gene
Will return a dictionary where each gene is a key and the value associated with each gene will be the pathway identifier

@param pathways, the dictionary of pathways
@return a dictionary with genes as keys and pathway identifiers as values
"""
def assignPathwayToGenes(pathways):
    allGenes = {}
    for pathwayIdentifier in pathways:
        pathway = pathways[pathwayIdentifier]
        for gene in pathway:
            allGenes[gene.strip()] = pathwayIdentifier
    return allGenes

"""
Compare each gene pair to the pathway data.

@param clusters, the dictionary of clusters
@param pathways, the dictionary of pathways
@return a confusion matrix as an array, [true positives, true negatives, false positives, false negatives]
"""
def compareClustersToPathways(clusters, pathways):
    confusion = [0, 0, 0, 0]
    allGenesFromClusters, clusterIdentifiers = assignClusterToGenes(clusters)
    allGenesFromPathways = assignPathwayToGenes(pathways)
    firstGeneIndex = 0
    secondGeneIndex = 1
    totalGenePairs = 0
    genePairsNotInPathwayData = 0
    print "Number of Genes From Clusters:", len(allGenesFromClusters)
    while firstGeneIndex < len(allGenesFromClusters):
        while secondGeneIndex < len(allGenesFromClusters):
            totalGenePairs+=1
            sameCluster = False
            samePathway = False
            if clusterIdentifiers[firstGeneIndex] == clusterIdentifiers[secondGeneIndex]:
                sameCluster = True
            gene1 = allGenesFromClusters[firstGeneIndex]
            gene2 = allGenesFromClusters[secondGeneIndex]
            #if either of the genes in the pair are not in the pathway data
            if gene1 not in allGenesFromPathways:
                secondGeneIndex+=1
                genePairsNotInPathwayData+=1
                continue
            if gene2 not in allGenesFromPathways:
                secondGeneIndex+=1
                genePairsNotInPathwayData+=1
                continue
            if allGenesFromPathways[gene1] == allGenesFromPathways[gene2]:
                samePathway = True
            if sameCluster == samePathway:
                if sameCluster:
                    confusion[0]+=1
                else:
                    confusion[1]+=1
            else:
                if sameCluster:
                    confusion[2]+=1
                else:
                    confusion[3]+=1
            secondGeneIndex+=1
        firstGeneIndex+=1
        secondGeneIndex = firstGeneIndex + 1
    print "Number of Genes Pairs Not In Pathway Data:", genePairsNotInPathwayData
    print "Number of Total Gene Pairs:", totalGenePairs
    print "Ratio of Gene Pairs Removed to Total Gene Pairs:", float(genePairsNotInPathwayData)/float(totalGenePairs)
    return confusion

"""
Calculates the analytic statistics.

@param confusion, the confusion matrix as an array, [true positives, true negatives, false positives, false negatives]
@return the accuracy, sensitivity, specificity, positive predictive value, and negative predictive value calculated from the confusion matrix
"""
def calcAnalyticStats(confusion):
    tp = confusion[0]
    tn = confusion[1]
    fp = confusion[2]
    fn = confusion[3]
    p = tp + fn
    n = fp + tn
    try:
        acc = float((tp+tn))/float((p+n))
    except ZeroDivisionError:
        acc = 0
    try:
        tpr = tp/float(p)
    except ZeroDivisionError:
        tpr = 0
    try:
        tnr = tn/float(n)
    except ZeroDivisionError:
        tnr = 0
    try:
        ppv = tp/float((tp+fp))
    except ZeroDivisionError:
        ppv = 0
    try:
        npv = tn/float((tn+fn))
    except ZeroDivisionError:
        npv = 0
    return acc, tpr, tnr, ppv, npv

"""
Display some analytic statistics about the comparison between the clusters and pathways
given the confusion matrix resulting from the comparison

@param confusion, the confusion matrix to calculate statistics from
"""
def displayAnalyticStats(confusion):
    tp = confusion[0]
    tn = confusion[1]
    fp = confusion[2]
    fn = confusion[3]
    acc, tpr, tnr, ppv, npv = calcAnalyticStats(confusion)
    print str(tp) + "\t|\t" + str(fp)
    print "---------------------"
    print str(fn) + "\t|\t" + str(tn)
    print "Accuracy:", acc
    print "Sensitivity:", tpr
    print "Specificity:", tnr
    print "Positive Predictive Value:", ppv
    print "Negative Predictive Value:", npv

"""
Appends the run times to a file for efficiency analysis

@param sampleSize, the size of the sample for this run instance
@param graphTime, the amount of time it takes to create the graph
@param clusterTime, the amount of time it takes to cluster the graph
@param compareTime, the amount of time it takes to compare clusters to pathways
"""
def appendTimesToFile(sampleSize, graphTime, clusterTime, compareTime):
    file = open("times.csv", "a")
    file.write(str(sampleSize)+", "+str(graphTime)+", "+str(clusterTime)+", "+str(compareTime)+"\n")
    file.close()

"""
Appends the analytic statistics to a file for analysis
@param sampleSize, the size of the sample for this run instance
@param confusion, the confusion matrix containing output information
"""
def appendDataToFile(sampleSize, confusion):
    file = open("data.csv", "a")
    acc, tpr, tnr, ppv, npv = calcAnalyticStats(confusion)
    file.write(str(sampleSize)+", "+str(acc)+", "+str(tpr)+", "+str(tnr)+", "+str(ppv)+", "+str(npv)+"\n")
    file.close()

"""
for i in range(100, 1000):
    begin = time.time()
    graph = sampleData("interactions.tsv", i)
    end = time.time()
    graphTime = end - begin
    begin = time.time()
    M, clusters = networkx_mcl(graph, inflate_factor = 1.8, max_loop = 60)
    clusters = convertClustersOfNodeIDsToClustersOfGeneIdentifiers(graph, clusters)
    clusters = removeDuplicateClusters(clusters)
    end = time.time()
    clusterTime = end - begin
    begin = time.time()
    pathways = readPathwayData("pathways.tsv")
    confusion = compareClustersToPathways(clusters, pathways)
    appendDataToFile(i, confusion)
    end = time.time()
    compareTime = end - begin
    appendTimesToFile(i, graphTime, clusterTime, compareTime)
    print i
"""

print "----------------------------------------------------------"

begin = time.time()
#graph = readDataToGraph("interactions.tsv")
originalGraph = sampleData("interactions.tsv", 1000)
graph = originalGraph
print "ORIGINAL GRAPH:"
#print "finished extracting interactions"
displayGraphStats(graph)
#print "displayed graph stats"
end = time.time()
print "Graph Creation Time:", end - begin

print ""

begin = time.time()
#nodes in cluster are identified by the same number as they were given
M, clusters = networkx_mcl(graph, inflate_factor = 1.8, max_loop = 60)
#print "finished clustering algorithm"
clusters = convertClustersOfNodeIDsToClustersOfGeneIdentifiers(graph, clusters)
clusters = removeDuplicateClusters(clusters)
displayClusterStats(clusters)
#print "displayed cluster stats"
end = time.time()
print "Clustering Time:", end - begin

print ""

begin = time.time()
pathways = readPathwayData("pathways.tsv")
#print "finished extracting pathways"
displayPathwayStats(pathways)
#print "displayed pathway stats"
confusion = compareClustersToPathways(clusters, pathways)
#print "finished validating clusters"

print ""

displayAnalyticStats(confusion)
#print "displayed analytic stats"
end = time.time()
print "Analysis Time:", end - begin

print "----------------------------------------------------------"

begin = time.time()
graph = generateRandomGraph(originalGraph)
print "RANDOM GRAPH:"
#print "finished extracting interactions"
displayGraphStats(graph)
#print "displayed graph stats"
end = time.time()
print "Graph Creation Time:", end - begin

print ""

begin = time.time()
#nodes in cluster are identified by the same number as they were given
M, clusters = networkx_mcl(graph, inflate_factor = 1.8, max_loop = 60)
#print "finished clustering algorithm"
clusters = convertClustersOfNodeIDsToClustersOfGeneIdentifiers(graph, clusters)
clusters = removeDuplicateClusters(clusters)
displayClusterStats(clusters)
#print "displayed cluster stats"
end = time.time()
print "Clustering Time:", end - begin

print ""

begin = time.time()
pathways = readPathwayData("pathways.tsv")
#print "finished extracting pathways"
displayPathwayStats(pathways)
#print "displayed pathway stats"
confusion = compareClustersToPathways(clusters, pathways)
#print "finished validating clusters"

print ""

displayAnalyticStats(confusion)
#print "displayed analytic stats"
end = time.time()
print "Analysis Time:", end - begin

print "----------------------------------------------------------"

begin = time.time()
graph = generateCompleteGraph(originalGraph)
print "COMPLETE GRAPH:"
#print "finished extracting interactions"
displayGraphStats(graph)
#print "displayed graph stats"
end = time.time()
print "Graph Creation Time:", end - begin

print ""

begin = time.time()
#nodes in cluster are identified by the same number as they were given
M, clusters = networkx_mcl(graph, inflate_factor = 1.8, max_loop = 60)
#print "finished clustering algorithm"
clusters = convertClustersOfNodeIDsToClustersOfGeneIdentifiers(graph, clusters)
clusters = removeDuplicateClusters(clusters)
displayClusterStats(clusters)
#print "displayed cluster stats"
end = time.time()
print "Clustering Time:", end - begin

print ""

begin = time.time()
pathways = readPathwayData("pathways.tsv")
#print "finished extracting pathways"
displayPathwayStats(pathways)
#print "displayed pathway stats"
confusion = compareClustersToPathways(clusters, pathways)
#print "finished validating clusters"

print ""

displayAnalyticStats(confusion)
#print "displayed analytic stats"
end = time.time()
print "Analysis Time:", end - begin

print "----------------------------------------------------------"

#nx.draw_circular(graph)
#plt.show()

#    clusters = dict with keys = [<cluster id>] values = [<vertex id>]
#    M = output matrix

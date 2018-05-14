import networkx as nx
import matplotlib.pyplot as plt
from random import shuffle

filename = "facebook_combined.txt"
#filename = "test.txt"
with open(filename) as file:
    DS = [tuple(line.rstrip('\n').split(' ')) for line in file]
shuffle(DS)
DSLen = len(DS)
DSDiv = int(DSLen*0.9)
G = nx.Graph()
G.add_edges_from(DS[0:DSDiv])
GTest = nx.Graph()
GTest.add_edges_from(DS[DSDiv:DSLen])
print("Datasets created")

##Adamic-Adar------------------------------------------------------------
GTestScore = nx.Graph()
for i in nx.adamic_adar_index(G):
    if i[2] > 0:
        if GTest.has_edge(i[0],i[1]) or GTest.has_edge(i[1],i[0]):
            GTestScore.add_edge(i[0],i[1],score=i[2])
print("Precision: ",GTestScore.number_of_edges()/GTest.number_of_edges())

##Common Neighbors--------------------------------------------------------
#GScore = nx.Graph()
#GTestScore = nx.Graph()
#nodeList = G.nodes()
#for n in nodeList:
#    nonNeighbors = nx.non_neighbors(G,n)
#    for nn in nonNeighbors:
#        scoreCN = len(sorted(nx.common_neighbors(G,n,nn)))
#        if scoreCN > 0:
#            GScore.add_edge(n,nn,score=scoreCN)
#            if (GTest.has_edge(n,nn) or GTest.has_edge(nn,n)) and not GTestScore.has_edge(n,nn) and not GTestScore.has_edge(nn,n):
#                GTestScore.add_edge(n,nn,score=scoreCN)
#print("Precision: ",GTestScore.number_of_edges()/GTest.number_of_edges())


#nx.draw(G, with_labels=True, font_weight='bold')
#plt.show()

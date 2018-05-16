#import matplotlib.pyplot as plt
import networkx as nx
from random import shuffle
import time
import itertools
import numpy
import scipy.linalg
import copy
import sys
from collections import defaultdict

def alg_katz(G,GTest):
    startTime = time.time() #BEGIN
    s,n = katz(G)
    print(s,n)
    finishTime = time.time() - startTime #END

    #print("Precision: ",GTestScore.number_of_edges()/GTest.number_of_edges())
    #print("Tiempo: ",time.time() - start_time)

def alg_simrank(G,GTest):
    startTime = time.time() #BEGIN
    s = simrank(G)
    print(s)
    finishTime = time.time() - startTime #END

    #print("Precision: ",GTestScore.number_of_edges()/GTest.number_of_edges())
    #print("Tiempo: ",time.time() - start_time)


def alg_adamic_adar(G,GTest):
    startTime = time.time() #BEGIN
    GTestScore = nx.Graph()
    for i in nx.adamic_adar_index(G):
        if i[2] > 0:
            if GTest.has_edge(i[0],i[1]):
                GTestScoreAA.add_edge(i[0],i[1],score=i[2])
    finishTime = time.time() - startTime #END

    print("Precision: ",GTestScore.number_of_edges()/GTest.number_of_edges())
    print("Tiempo: ",finishTime)

def alg_common_neighbors(G,GTest):
    startTime = time.time() #BEGIN
    GScore = nx.Graph()
    GTestScore = nx.Graph()
    nodeList = sorted(G.nodes())
    for n in nodeList:
        for nn in nx.non_neighbors(G,n):
            scoreCN = len(sorted(nx.common_neighbors(G,n,nn)))
            if 0 < scoreCN:
                GScore.add_edge(n,nn,score=scoreCN)
                if GTest.has_edge(n,nn):
                    GTestScore.add_edge(n,nn,score=scoreCN)
    finishTime = time.time() - startTime #END

    print("Precision: ",GTestScore.number_of_edges()/GTest.number_of_edges())
    print("Tiempo: ",finishTime)

##-------------------------------------------------------------------------------
def katz(G, c=0.9, remove_neighbors=False, inv_method=0):
  if type(G) == nx.MultiGraph or type(G) == nx.MultiDiGraph:
    raise Exception("katz() not defined for graphs with multiedges.")

  if G.is_directed():
    raise Exception("katz() not defined for directed graphs.")

  A = nx.adjacency_matrix(G, nodelist=G.nodes(), weight=None)
  w, v = numpy.linalg.eigh(A.toarray())
  lambda1 = max([abs(x) for x in w])
  I = numpy.eye(A.shape[0])
  S = None
  if inv_method == 1:
    S = scipy.linalg.pinv(I - c/lambda1 * A)
  elif inv_method == 2:
    S = numpy.linalg.inv(I - c/lambda1 * A)
  else:
    S = numpy.linalg.pinv(I - c/lambda1 * A)
  return S, G.nodes()


def simrank(G, c=0.9, max_iter=100, remove_neighbors=False, remove_self=False, dump_process=False):
  if type(G) == nx.MultiGraph or type(G) == nx.MultiDiGraph:
    raise Exception("simrank() not defined for graphs with multiedges.")

  if G.is_directed():
    raise Exception("simrank() not defined for directed graphs.")

  sim_old = defaultdict(list)
  sim = defaultdict(list)
  for n in G.nodes():
    sim[n] = defaultdict(int)
    sim[n][n] = 1
    sim_old[n] = defaultdict(int)
    sim_old[n][n] = 0

  # calculate simrank
  for iter_ctr in range(max_iter):
    if _is_converge(sim, sim_old):
      break
    sim_old = copy.deepcopy(sim)
    for i, u in enumerate(G.nodes()):
      if dump_process:
        sys.stdout.write("\r%d : % d / %d" % (iter_ctr, i, G.number_of_nodes()))
      for v in G.nodes():
        if u == v:
          continue
        s_uv = 0.0
        for n_u in G.neighbors(u):
          for n_v in G.neighbors(v):
            s_uv += sim_old[n_u][n_v]
        sim[u][v] = (c * s_uv / (len(list(G.neighbors(u))) * len(list(G.neighbors(v))))) \
            if len(list(G.neighbors(u))) * len(list(G.neighbors(v))) > 0 else 0
    if dump_process:
      print('')

  if remove_self:
    for m in G.nodes():
      G[m][m] = 0

  if remove_neighbors:
    for m in G.nodes():
      for n in G.neighbors(m):
        sim[m][n] = 0

  return sim

def _is_converge(s1, s2, eps=1e-4):
  for i in s1.keys():
    for j in s1[i].keys():
      if abs(s1[i][j] - s2[i][j]) >= eps:
        return False
  return True
##------------------------------------------------------------------------------

if __name__ == "__main__":
    #filename = "facebook_combined.txt"
    filename = "test.txt"
    with open(filename) as file:
        DS = [tuple(line.rstrip('\n').split(' ')) for line in file]
    shuffle(DS)
    DSLen = len(DS)
    DSDiv = int(DSLen*0.9)
    G = nx.Graph()
    G.add_edges_from(DS[0:DSDiv])
    GTest = nx.Graph()
    GTest.add_edges_from(DS[DSDiv:DSLen])
    alg_simrank(G,GTest)

#nx.draw(G, with_labels=True, font_weight='bold')
#plt.show()

#import matplotlib.pyplot as plt
import networkx as nx
from random import shuffle, randint
import time
import itertools
import numpy
import scipy.linalg
import copy
import sys
from collections import defaultdict

ID = randint(10000,99999)

def alg_katz(G,GTest):
    startTime = time.time() #BEGIN
    GScoreRight = nx.Graph()
    GScoreWrong = nx.Graph()
    v = list(G.nodes())
    nv = len(v)
    k = katz(G)
    for i in range(0,nv):
        for j in range(0,nv):
            vi = v[i]
            vj = v[j]
            scoreK = k.item((i,j))
            if scoreK > 0.000001 and not i == j and not G.has_edge(vi,vj) and not GScoreRight.has_edge(vi,vj) and not GScoreWrong.has_edge(vi,vj):
                if GTest.has_edge(vi,vj):
                    GScoreRight.add_edge(vi,vj,score=scoreK)
                else:
                    GScoreWrong.add_edge(vi,vj,score=scoreK)
    finishTime = time.time() - startTime #END
    measures("KATZ", GTest, GScoreRight, GScoreWrong, finishTime)

def alg_simrank(G,GTest):
    startTime = time.time() #BEGIN
    GScoreRight = nx.Graph()
    GScoreWrong = nx.Graph()
    for t1,t2 in simrank(G).items():
        for sc in t2.items():
            sc0 = sc[0]
            scoreSR = sc[1]
            if scoreSR > 0.000001 and not t1 == sc0 and not G.has_edge(t1,sc0) and not GScoreRight.has_edge(t1,sc0) and not GScoreWrong.has_edge(t1,sc0):
                if GTest.has_edge(t1,sc0):
                    GScoreRight.add_edge(t1,sc0,score=scoreSR)
                else:
                    GScoreWrong.add_edge(t1,sc0,score=scoreSR)
    finishTime = time.time() - startTime #END
    measures("SIMRANK", GTest, GScoreRight, GScoreWrong, finishTime)

def alg_adamic_adar(G,GTest):
    startTime = time.time() #BEGIN
    GScoreRight = nx.Graph()
    GScoreWrong = nx.Graph()
    for i in nx.adamic_adar_index(G):
        i0 = i[0]
        i1 = i[1]
        i2 = i[2]
        if i[2] > 0:
            if GTest.has_edge(i0,i1):
                GScoreRight.add_edge(i0,i1,score=i2)
            else:
                GScoreWrong.add_edge(i0,i1,score=i2)
    finishTime = time.time() - startTime #END
    measures("ADAMIC-ADAR", GTest, GScoreRight, GScoreWrong, finishTime)

def alg_common_neighbors(G,GTest):
    startTime = time.time() #BEGIN
    GScoreRight = nx.Graph()
    GScoreWrong = nx.Graph()
    nodeList = list(G.nodes())
    for n in nodeList:
        for nn in nx.non_neighbors(G,n):
            scoreCN = len(sorted(nx.common_neighbors(G,n,nn)))
            if 0 < scoreCN  and not GScoreRight.has_edge(n,nn) and not GScoreWrong.has_edge(n,nn):
                if GTest.has_edge(n,nn):
                    GScoreRight.add_edge(n,nn,score=scoreCN)
                else:
                    GScoreWrong.add_edge(n,nn,score=scoreCN)
    finishTime = time.time() - startTime #END
    measures("COMMON NEIGHBORS", GTest, GScoreRight, GScoreWrong, finishTime)

def measures(ALG, GTest, GScoreRight, GScoreWrong, finishTime):
    PRECISION = GScoreRight.number_of_edges()/GTest.number_of_edges()
    TIEMPO = finishTime
    
    n = int(GScoreRight.number_of_edges() * 0.9)
    ni = 0
    nii = 0
    lisR = list(GScoreRight.edges())
    lis_r = list(range(len(lisR)))
    shuffle(lis_r)
    lisW = list(GScoreWrong.edges())
    lis_w = list(range(len(lisW)))
    shuffle(lis_w)
    for i in range(0,n):
        if GScoreRight.get_edge_data(*lisR[lis_r[i]])['score'] > GScoreWrong.get_edge_data(*lisW[lis_w[i]])['score']:
            ni += 1
        else:
            nii += 1
    AUC = (ni + (nii * 0.5)) / n
    
    file = open(str(ID)+"_RESULT"+".txt","a")
    file.write(ALG+":")
    file.write("\n")
    file.write("  PRECISION: "+str(PRECISION)) 
    file.write("\n")
    file.write("  TIEMPO: "+str(TIEMPO)) 
    file.write("\n")
    file.write("  AUC: "+str(AUC))
    file.write("\n")
    file.write("\n")
    file.close()

    nx.write_edgelist(GScoreRight,str(ID)+"_GR_"+ALG,data=['score'])
    nx.write_edgelist(GScoreWrong,str(ID)+"_GW_"+ALG,data=['score'])

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
  return S


def simrank(G, c=0.9, max_iter=1, remove_neighbors=False, remove_self=False, dump_process=False):
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
    nx.write_edgelist(G,str(ID)+"_G")
    nx.write_edgelist(GTest,str(ID)+"_GTest")
    alg_common_neighbors(G,GTest)
    alg_adamic_adar(G,GTest)
    alg_katz(G,GTest)
    alg_simrank(G,GTest)

#nx.draw(G, with_labels=True, font_weight='bold')
#plt.show()

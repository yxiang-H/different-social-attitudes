import numpy as np
import random
import networkx as nx
from AD import PowerLaw_distribution as dis
class simplex_vertex:
    def __init__(self, vertex_id, activity):
        self.act = activity
        self.name = vertex_id

def add_k_edges(n, e, tgraph):  # 往图中添加边
    new_neigh = map(lambda x: [n, x], e)
    tgraph.add_edges_from(new_neigh)
    return tgraph


class model(object):
    def __init__(self, vertex_dict, m , eta, del_t, n, T, p, delta, beta, mu, last_t, I_0):
        self.vertex_dict = vertex_dict
        self.m = m           #产生m条边
        self.eta = eta
        self.del_t = del_t
        self.n = n            #节点数
        self.T = T            #演化的时间步
        self.p = p            #A类人的比例
        self.delta = delta    #减少的连边比例
        self.delta_m = int(delta * m)
        self.last_t = last_t  #最后记录的时间步
        self.I_0 = I_0        #初始感染种子数
        self.beta = beta      #感染率
        self.mu = mu          #恢复率
    def initial(self):
        sAgent = set()
        iAgent = set()
        asAgent = set()   #a 为风险者，不会减少连边
        bsAgent = set()   #b 为担忧者，会减少连边
        A_0 = int(self.p * self.n)   #A类人有p*n
        for node in self.vertex_dict.keys():
            sAgent.add(node)
        for _ in range(I_0):
            to_infect = random.choice(list(sAgent))
            iAgent.add(to_infect)
            sAgent.remove(to_infect)
        for bnodes in self.vertex_dict.keys():
            bsAgent.add(bnodes)
        for _ in range(A_0):
            to_a = random.choice(list(bsAgent))
            asAgent.add(to_a)
            bsAgent.remove(to_a)
        aiAgent = asAgent & iAgent
        biAgent = bsAgent & iAgent
        asAgent = asAgent - aiAgent
        bsAgent = bsAgent - biAgent
        return asAgent, bsAgent, aiAgent, biAgent
    def reduce_m(self):
        x = self.initial()
        forward_asnodes = x[0]
        forward_bsnodes = x[1]
        forward_ainodes = x[2]
        forward_binodes = x[3]
        Ai_frac = 0
        Bi_frac = 0
        for t in range(self.T):
            if not forward_ainodes|forward_binodes:
                break
            forward_as = set()
            forward_bs = set()
            forward_ai = set()
            forward_bi = set()
            tgraph = nx.Graph()
            tgraph.add_nodes_from(self.vertex_dict.keys())
            for node in forward_asnodes:
                tgraph.nodes[node]['type'] = 'A'
                tgraph.nodes[node]['state'] = 'S'
                tgraph.nodes[node]['activity'] = self.vertex_dict[node].act
                forward_as.add(node)
            for i in forward_bsnodes:
                tgraph.nodes[i]['type'] = 'B'
                tgraph.nodes[i]['state'] = 'S'
                tgraph.nodes[i]['activity'] = self.vertex_dict[i].act
                forward_bs.add(i)
            for a in forward_ainodes:
                tgraph.nodes[a]['type'] = 'A'
                tgraph.nodes[a]['state'] = 'I'
                tgraph.nodes[a]['activity'] = self.vertex_dict[a].act
                forward_ai.add(a)
            for b in forward_binodes:
                tgraph.nodes[b]['type'] = 'B'
                tgraph.nodes[b]['state'] = 'I'
                tgraph.nodes[b]['activity'] = self.vertex_dict[b].act
                forward_bi.add(b)
            for n1 in forward_asnodes|forward_ainodes:
                if np.random.rand() <= tgraph.nodes[n1]['activity'] * self.eta * self.del_t:
                    nodes = list(self.vertex_dict.keys())
                    nodes.remove(n1)
                    for i in sorted(nx.neighbors(tgraph, n1)):
                        nodes.remove(i)
                    neigh = random.sample(nodes, self.m)
                    tgraph = add_k_edges(n1, neigh, tgraph)
            for n2 in forward_bsnodes|forward_binodes:
                if np.random.rand() <= tgraph.nodes[n2]['activity'] * self.eta * self.del_t:
                    nodes = list(self.vertex_dict.keys())
                    nodes.remove(n2)
                    for i in sorted(nx.neighbors(tgraph, n2)):
                        nodes.remove(i)
                    neigh = random.sample(nodes, self.delta_m)
                    tgraph = add_k_edges(n2, neigh, tgraph)
            for node_as in forward_asnodes:
                node_s_neighbors = sorted(nx.neighbors(tgraph, node_as))
                i = 0
                infected = False
                while i < len(node_s_neighbors):
                    if tgraph.nodes[node_s_neighbors[i]]['state'] == 'I':
                        infected = np.random.random() <= self.beta
                        if infected == True: break
                    i += 1
                if infected:
                    forward_ai.add(node_as)
                    forward_as.remove(node_as)
            for node_bs in forward_bsnodes:
                node_s_neighbors = sorted(nx.neighbors(tgraph, node_bs))
                i = 0
                infected = False
                while i < len(node_s_neighbors):
                    if tgraph.nodes[node_s_neighbors[i]]['state'] == 'I':
                        infected = np.random.random() <= self.beta
                        if infected == True: break
                    i += 1
                if infected:
                    forward_bi.add(node_bs)
                    forward_bs.remove(node_bs)
            for node_ai in forward_ainodes:
                if np.random.random() <= self.mu:
                    forward_as.add(node_ai)
                    forward_ai.remove(node_ai)
            for node_bi in forward_binodes:
                if np.random.random() <= self.mu:
                    forward_bs.add(node_bi)
                    forward_bi.remove(node_bi)
            forward_asnodes = forward_as
            forward_bsnodes = forward_bs
            forward_ainodes = forward_ai
            forward_binodes = forward_bi
            if t >= self.last_t:
                ai_simulation = float(len(forward_ainodes)) / self.n
                Ai_frac += ai_simulation
                bi_simulation = float(len(forward_binodes)) / self.n
                Bi_frac += bi_simulation
        return Ai_frac/(self.T-self.last_t), Bi_frac/(self.T-self.last_t)


def pro_acts(N):
    dist = []
    for i in range(N):
        dist.append(dis.inv_cdf(2.8, 0.01)(random.uniform(0, 1)))
    return dist

beta_s = [0.02 * i for i in range(51)]
p_s = [0.02 * i for i in range(51)]
# repeat = 1
N = pow(10,4)
n = 5000
m = 8
eta = 1
del_t = 1
T = 10000
delta = 0.5
mu = 0.1
last_t = 9800
I_0 = 50

###重复多次实验
# for beta in beta_s:
#     stationary_AI = []
#     stationary_BI = []
#     for p in p_s:
#         for i in range(repeat):
#             dist = pro_acts(N)
#             act = np.random.choice(dist, n)
#             dist = []
#             vertex_dict = {}
#             for j in range(n):
#                 vertex_dict[j] = simplex_vertex(j, act[j])
#             uau_sis = model(vertex_dict, m=m, eta=eta, del_t=del_t, n=n, T=T, p=p, delta=delta, beta=beta, mu=mu, last_t=last_t, I_0=I_0)
#             stationary_AI.append(uau_sis.reduce_m()[0])
#             stationary_BI.append(uau_sis.reduce_m()[1])
#         mean_AI = sum(stationary_AI)/len(stationary_AI)
#         mean_BI = sum(stationary_BI)/len(stationary_BI)
#         with open('D:/output1/frequency.txt', 'a') as f:
#             f.write(str(beta)+' '+ str(p)+' '+ str(mean_AI)+' '+str(mean_BI)+'\n')


### 重复一次实验
for beta in beta_s:
    for p in p_s:
        dist = pro_acts(N)
        act = np.random.choice(dist, n)
        dist = []
        vertex_dict = {}
        for j in range(n):
            vertex_dict[j] = simplex_vertex(j, act[j])
        uau_sis = model(vertex_dict, m=m, eta=eta, del_t=del_t, n=n, T=T, p=p, delta=delta, beta=beta, mu=mu,
                        last_t=last_t, I_0=I_0)
        mean_AI, mean_BI = uau_sis.reduce_m()
        with open('D:/output1/frequency.txt', 'a') as f:
            f.write(str(beta)+' '+ str(p)+' '+ str(mean_AI)+' '+str(mean_BI)+'\n')













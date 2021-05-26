###############################################################################

##### --- Author:  Hend Alrasheed --- #####
##### --- Date:    May, 2021      --- #####
##### --- Purpose: This code analyzes the shortest paths utilization  --- #####
##### ---          in 16 real and artificial networks during the IM   --- #####
##### ---          (influence maximization) algorithm. Three seed     --- #####
##### ---          selection algorithms are used: random, degree,     --- #####
##### ---          and closeness. The IC model is used for diffusion. --- #####

###############################################################################


import numpy as np
import matplotlib.pyplot as plt
import random
import networkx as nx
import scipy
import matplotlib.ticker as mtick
import matplotlib.patches as mpatches
import math
from matplotlib.lines import Line2D
from matplotlib import gridspec

#-----------------------------------------------------------


def choose_seeds_random(G, num):
    
    seeds = random.sample(G.nodes(), num)
    for node in seeds:
        G.nodes[node]['state'] = 'A'
        G.nodes[node]['whichSeed'] = node
        G.nodes[node]['seed'] = 'true'

    #print(str(num) + " random nodes have been assigned as seed nodes...")
    #print 'num = ', len([n for n,d in G.nodes(data=True) if d['state']=='A'])

#=============================================

def choose_seeds_degree(G, num):

    s = sorted(G.degree, key=lambda x: x[1], reverse=True)
    seeds = s[:num]
    for key,value in seeds:
        G.nodes[key]['state'] = 'A'
        G.nodes[key]['whichSeed'] = key
        G.nodes[key]['seed'] = 'true'

#=============================================

def choose_seeds_closeness(G, num): 

    c = nx.closeness_centrality(G)
    sorted_c = sorted(c, key=c.get, reverse=True)
    print("\nTop closeness cent nodes:")
    seeds = list(sorted_c)
    print seeds[:num]
        
    for node in seeds:
        G.nodes[node]['state'] = 'A'
        G.nodes[node]['whichSeed'] = node
        G.nodes[node]['seed'] = 'true'

#============================================    

def activation(G):

    num=0
    seeds = [n for n,d in G.nodes(data=True) if d['state']=='A']
   
    for node in seeds:
            for neighbor in G.neighbors(node):
                if G.nodes[neighbor]['state'] == 'N':
                   if random.random() < p:
                      num=num+1 
                      G.nodes[neighbor]['state'] = 'A'
                      G.nodes[neighbor]['whichSeed'] = G.nodes[node]['whichSeed']
                      G.nodes[neighbor]['infCounter'] = G.nodes[node]['infCounter'] + 1
            G.nodes[node]['state'] = 'P'          

    return num

#=============================================                      

def shortest_path_analysis(G):

    total = shortest = 0
    diff = [0] * maxDiff
    
    for node in G.nodes():
        if G.nodes[node]['state'] == 'A' or G.nodes[node]['state'] == 'P':
            total+=1
            d = nx.shortest_path_length(G, node, G.nodes[node]['whichSeed'])
            if d == G.nodes[node]['infCounter']:
                shortest+=1
            diff[G.nodes[node]['infCounter']-d]+=1
 
    return diff

#=============================================

def shortest_path_analysis_by_distance(G):

    total = shortest = 0
    diff = [0] * maxDiff
    
    for node in G.nodes():
        if G.nodes[node]['state'] == 'A' or G.nodes[node]['state'] == 'P':
            d = nx.shortest_path_length(G, node, G.nodes[node]['whichSeed'])
            if d != G.nodes[node]['infCounter']:
               diff[d]+=1

    return diff

#=============================================

def shortest_path_analysis_by_eccentricity(G):

    radius = nx.radius(G)
    total = shortest = 0
    diff = [0] * maxDiff
    
    for node in G.nodes():
        if G.nodes[node]['state'] == 'A' or G.nodes[node]['state'] == 'P':
            d = nx.shortest_path_length(G, node, G.nodes[node]['whichSeed'])
            if d != G.nodes[node]['infCounter']:
               e = nx.eccentricity(G,v=node) 
               diff[e-radius]+=1

    return diff

#=============================================
   
def IM_simulation(G, i, seed_method, name):

    infection_history=[]
    susceptiblity_history=[]
    outbreak_size = 0

    nx.set_node_attributes(G, 'N', name='state') # node state can be: Active (A), Non-active (N), or Processed (P). Active nodes represent nodes who are currently active. 
    nx.set_node_attributes(G, 'false', name='seed')
    nx.set_node_attributes(G, -1000, name='whichSeed')
    nx.set_node_attributes(G, 0, name='infCounter')

    num_of_seeds = int(math.ceil(percent_of_seed_nodes*len(G)))
    
    if seed_method == 1:
        choose_seeds_random(G,num_of_seeds)
    elif seed_method == 2:
        choose_seeds_degree(G,num_of_seeds)
    elif seed_method == 3:
        choose_seeds_closeness(G,num_of_seeds)

    infection_history.append(num_of_seeds)

    for i in range(num_of_time_steps-1):
        new_infections = activation(G)
        infection_history.append(infection_history[-1] + new_infections)
        outbreak_size += new_infections

    differences = shortest_path_analysis(G)
    differencesByDist = shortest_path_analysis_by_distance(G)
    differencesByEcc = shortest_path_analysis_by_eccentricity(G)
    
    return differences, differencesByDist, differencesByEcc, outbreak_size

    
#=============================================

def IM(G, name, seed_method):

    global unutlized_list
    results = [0] * maxDiff
    resultsDist = [0] * maxDiff
    resultsEcc = [0] * maxDiff
    OB = 0

    for i in range(num_of_simulations):
        A, B, C, OS = IM_simulation(G,i,seed_method, name)
        results = [x + y for x, y in zip(results, A)]
        resultsDist = [x + y for x, y in zip(resultsDist, B)]
        resultsEcc = [x + y for x, y in zip(resultsEcc, C)]
        OS = len([n for n,d in G.nodes(data=True) if d['state']=='A' or d['state']=='P'])
        OB += OS
        

    results = [round(float(i)/num_of_simulations) for i in results]
    resultsDist = [round(float(i)/num_of_simulations) for i in resultsDist]
    resultsEcc = [round(float(i)/num_of_simulations) for i in resultsEcc]
    outbreak_size = round(OB/num_of_simulations)

    total=0
    for i in range(1,len(results)):
        total += results[i]

    totalD=0
    for i in range(1,len(resultsDist)):
        totalD += resultsDist[i]

    totalE=0
    for i in range(1,len(resultsEcc)):
        totalE += resultsEcc[i]

    if seed_method == 1:
        t='Random'
    elif seed_method == 2:
        t='Degree'
    elif seed_method == 3:
        t='Closeness'
    

    print 'Graph: ', name
    print 'Method: ', t
    print 'Avg outbreak size: ', outbreak_size/len(G)*100
    print '% of un-utilized shortest paths = ', total/(total+results[0])*100    
    unutlized_list.append(total/(total+results[0])*100)
    return results, resultsDist, resultsEcc, outbreak_size/len(G)*100, total/(total+results[0])*100
    

#=============================================


graphs=[]


G1 = nx.karate_club_graph()
G2 = nx.read_edgelist('dataset/emailcc.txt')
G3 = nx.read_edgelist('dataset/facebook.txt')
G4 = nx.read_edgelist('dataset/DutchEliteCC.txt')
G5 = nx.read_edgelist('dataset/wattsStrogatz.3.txt')
G6 = nx.read_edgelist('dataset/wattsStrogatz.5.txt')
G7 = nx.read_edgelist('dataset/wattsStrogatz.7.txt')
G8 = nx.read_edgelist('dataset/BarabasiAlbert3.txt')
G9 = nx.read_edgelist('dataset/BarabasiAlbert5.txt')
G10 = nx.read_edgelist('dataset/BarabasiAlbert7.txt')
G11 = nx.read_edgelist('dataset/ER_1582.txt')
G12 = nx.read_edgelist('dataset/ER_1993.txt')
G13 = nx.read_edgelist('dataset/ER_2498_0.0032.txt')
G14 = nx.read_edgelist('dataset/PL_2236-h4.txt')
G15 = nx.read_edgelist('dataset/PL_1761-h7.txt')
G16 = nx.read_edgelist('dataset/PL_1199-h9.txt')


graphs.append(G1)
graphs.append(G2)
graphs.append(G3)
graphs.append(G4)
graphs.append(G5)
graphs.append(G6)
graphs.append(G7)
graphs.append(G8)
graphs.append(G9)
graphs.append(G10)
graphs.append(G11)
graphs.append(G12)
graphs.append(G13)
graphs.append(G14)
graphs.append(G15)
graphs.append(G16)


maxDiff = 50
p = .25
num_of_simulations = 100
percent_of_seed_nodes = 0.2
num_of_time_steps = 100
unutlized_list=[]
names=["Karate \n Club","Email","Facebook","Dutch \n Elite","Watts \n Strogatz(0.3)","Watts \n Strogatz(0.5)","Watts \n Strogatz(0.7)","Barabasi \n Albert(3)",
"Barabasi \n Albert(5)","Barabasi \n Albert(7)","Erdos \n Renyi(1.6)","Erdos \n Renyi(2)","Erdos \n Renyi(8)","Power \n Law(1.8)","Power \n Law(2)","Power \n Law(2.7)"]
 
print('Start')

for i in range(len(graphs)):

    G = graphs[i]
    name = names[i]
    diffR, distR, eccR, obsizeR, unR = IM(G, name, 1)
    diffD, distD, eccD, obsizeD, unD = IM(G, name, 2)
    diffG, distG, eccG, obsizeG, unG = IM(G, name, 3)
    
    print '\n==============================================\n'
  
print '% of unutilized paths in each graph:'
print unutlized_list
print('Done!')


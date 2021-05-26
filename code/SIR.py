###############################################################################

##### --- Author:  Hend Alrasheed --- #####
##### --- Date:    May, 2021      --- #####
##### --- Purpose: This code analyzes the shortest paths utilization  --- #####
##### ---          in 16 real and artificial networks during the SIR  --- #####
##### ---          (susceptible-infectious-recovered) model.          --- #####

###############################################################################

import numpy as np
import matplotlib.pyplot as plt
import random
import networkx as nx
import scipy

#-----------------------------------------------------------



#=============================================

def choose_seeds(num,G):

    seeds = random.sample(G.nodes(),num)
    for node in seeds:
        G.nodes[node]['state'] = 'I'
        G.nodes[node]['whichSeed'] = node
        G.nodes[node]['seed'] = 'true'

#=============================================

def diffusion(time_step, p,G):

    num=0

    for node in G.nodes():
        if G.nodes[node]['state'] == 'I' and G.nodes[node]['infection_time'] != time_step:
            for neighbor in G.neighbors(node):
                if G.nodes[neighbor]['state'] == 'S':
                   if random.random() < p:
                      num=num+1 
                      G.nodes[neighbor]['state'] = 'I'
                      G.nodes[neighbor]['infection_time'] = time_step
                      G.nodes[neighbor]['whichSeed'] = G.nodes[node]['whichSeed']
                      G.nodes[neighbor]['infCounter'] = G.nodes[node]['infCounter'] + 1

    return num                       

#=============================================                      

def recovery(time_step, p, T,G):

    num=0
    
    for node in G.nodes():
        if G.nodes[node]['state'] == 'I' and ((time_step - G.nodes[node]['infection_time']) >= T):
           if random.random() < p:
              num=num+1
              G.nodes[node]['state'] = 'R'

    return num         


#=============================================

def shortest_path_analysis(G):

    total = shortest = 0
    diff = [0] * maxDiff
    
    for node in G.nodes():
        if G.nodes[node]['state'] == 'I' or G.nodes[node]['state'] == 'R':
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
        if G.nodes[node]['state'] == 'I' or G.nodes[node]['state'] == 'R':
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
        if G.nodes[node]['state'] == 'I' or G.nodes[node]['state'] == 'R':
            d = nx.shortest_path_length(G, node, G.nodes[node]['whichSeed'])
            if d != G.nodes[node]['infCounter']:
               e = nx.eccentricity(G,v=node) 
               diff[e-radius]+=1

    return diff

#=============================================

def generate_lists(diff):

    #remove zero elements
    index = 0
    for i in range(len(diff)-1):
        if diff[i]!=0 and diff[i+1]==0:
            index=i+1
    diff = diff[0:index]
    values = np.arange(0,len(diff))
    total=sum(diff)
    percents = [(x / total)*100 for x in diff]
    return values,diff,percents


#=============================================
    
def SIR_simulation(G, i):

    global infection_history, recovery_history, susceptiblity_history, outbreak_size

    nx.set_node_attributes(G, 'S', name='state') # node state can be: Susceptible (S), Infected (I), or Recovered (R)
    nx.set_node_attributes(G, int(0), name='infection_time')
    nx.set_node_attributes(G, 'false', name='seed')
    nx.set_node_attributes(G, -1000, name='whichSeed')
    nx.set_node_attributes(G, 0, name='infCounter')
    nx.set_node_attributes(G, 0, name='visitCounter')
    infection_history.append(num_of_seed_nodes)
    recovery_history.append(0)
    susceptiblity_history.append(len(G) - num_of_seed_nodes)
    choose_seeds(num_of_seed_nodes,G)

    for i in range(num_of_time_steps-1): 
        new_infections = diffusion(i, infection_rate,G)
        new_recoveries = recovery(i, recovery_rate, recovery_threshold,G) 
        infection_history.append(infection_history[-1] + new_infections)
        recovery_history.append(recovery_history[-1] + new_recoveries)
        susceptiblity_history.append(len(G) - sum(1 for n,d in G.nodes(data=True) if d['state']=='I') - sum(1 for n,d in G.nodes(data=True) if d['state']=='R'))
        outbreak_size += new_infections

    differences = shortest_path_analysis(G)
    differencesByDist = shortest_path_analysis_by_distance(G)
    differencesByEcc = shortest_path_analysis_by_eccentricity(G)
    
    return differences, differencesByDist, differencesByEcc

    
#=============================================

def SIR(G, name):

    global outbreak_size
    global unutlized_list
    results = [0] * maxDiff
    resultsDist = [0] * maxDiff
    resultsEcc = [0] * maxDiff

    for i in range(num_of_simulations):
        A, B, C = SIR_simulation(G,i)
        results = [x + y for x, y in zip(results, A)]
        resultsDist = [x + y for x, y in zip(resultsDist, B)]
        resultsEcc = [x + y for x, y in zip(resultsEcc, C)]

    results = [round(float(i)/num_of_simulations) for i in results]
    resultsDist = [round(float(i)/num_of_simulations) for i in resultsDist]
    resultsEcc = [round(float(i)/num_of_simulations) for i in resultsEcc]
    outbreak_size = round(outbreak_size/num_of_simulations)

    total=0
    for i in range(1,len(results)):
        total += results[i]

    total=0
    for i in range(1,len(resultsDist)):
        total += resultsDist[i]

    total=0
    for i in range(1,len(resultsEcc)):
        total += resultsEcc[i]     

    print 'Differences list: ', results
    print 'Graph: ', name    
    print 'Avg outbreak size: ', outbreak_size
    print '% of un-utilized shortest paths = ', total/(total+results[0])*100
    print '-------------------------------------'     
    unutlized_list.append(total/(total+results[0])*100)
    return results, resultsDist, resultsEcc, outbreak_size
    



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


maxDiff = 50
infection_rate=.5
recovery_rate=.1
recovery_threshold = 2
infection_history=[]
susceptiblity_history=[]
recovery_history=[]
num_of_simulations = 100
num_of_seed_nodes = 1
num_of_time_steps = 100
outbreak_size=0.0
G=nx.Graph()
unutlized_list=[]
print('Start')


results1, results21, results31, obsize = SIR(G1, "Karate-Club")
results2, results22, results32, obsize = SIR(G2, "Email")
results3, results23, results33, obsize = SIR(G3, "Facebook")
results4, results24, results34, obsize = SIR(G4, "Dutch-Elite")
results5, results25, results35, obsize = SIR(G5, "Watts-Strogatz(0.3)")
results6, results26, results36, obsize = SIR(G6, "Watts-Strogatz(0.5)")
results7, results27, results37, obsize = SIR(G7, "Watts-Strogatz(0.7)")
results8, results28, results38, obsize = SIR(G8, "Barabasi-Albert(3)")
results9, results29, results39, obsize = SIR(G9, "Barabasi-Albert(5)")
results10, results210, results310, obsize = SIR(G10, "Barabasi-Albert(7)")
results11, results211, results311, obsize = SIR(G11, "Erdos-Renyi(1.6)")
results12, results212, results312, obsize = SIR(G12, "Erdos-Renyi(2)")
results13, results213, results313, obsize = SIR(G13, "Erdos-Renyi(8)")
results14, results214, results314, obsize = SIR(G14, "Power-Law(1.8)")
results15, results215, results315, obsize = SIR(G15, "Power-Law(2)")
results16, results216, results316, obsize = SIR(G16, "Power-Law(2.7)")


print '% of unutilized paths in each graph:'
print unutlized_list

print('Done!')

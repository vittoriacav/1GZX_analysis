#!/usr/bin/env python3
import sys
import networkx as nx
import copy
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress

def parse_mitab(filename, org, id1 = 0, id2 = 1):
	gi = {}
	f = open (filename, 'r')
	for line in f:
		if line[0] == '#':
			continue
		line = line.rstrip()
		v = line.split('\t') 
		
		# Select only the lines in which both the interacting proteins 
		#  have a uniprotkb identifier
		# v[id1] = v[0] = unique identifier of interactor A (first column)
		# v[id2] = v[1] = unique identifier of interactor B (second column)
		
		if v[id1].find('uniprotkb:')==-1 or v[id2].find('uniprotkb:')==-1: 
			continue
		
		# We want consider only proteins interacting in a specific organism
		#  that in our case is human (taxid: 9606), hence we check the taxid
		#  for the 2 current interacting proteins
		# v[9] = NCBI taxonomy identifier for interactor A
		# v[10] = NCBI taxonomy identifier for interactor B
		
		if v[9].find('taxid:'+org)==-1 or v[10].find('taxid:'+org)==-1:
			continue
		
		# We split to obtain protein A and protein B names
		gid1 = v[id1].split('|')[0].split(':')[1].split('-')[0]
		gid2 = v[id2].split('|')[0].split(':')[1].split('-')[0]
		gids = [gid1, gid2]
		gids.sort()
		gi[tuple(gids)] = True
	return list(gi.keys())
	
if __name__ == '__main__':
	# As input we need to provide the filename of the database of 
	#  protein-protein interaction. In our case we will use "direct_intact.txt"
	#  that contains only direct interactions between proteins
	#  As second argument for the input we need to provide the TAXID of 
	#  the desired organism. In our case TAXID human: 9606
	filename = sys.argv[1]
	org = sys.argv[2] 
	gis = parse_mitab(filename, org)
	
#NOW WE CREATE THE GENERAL GRAPH OF DIRECT PROTEIN-PROTEIN INTERACTION 
	g = nx.Graph()
	g.add_edges_from(gis)
	
	# To plot the whole graph based on gis
	#nx.draw(g)
	#plt.show()

#NOW WE CREATE THE GRAPH of subAlpha and subBeta

	j = 0
	max_len = 0
	sum_nodes = 0
	for i in nx.connected_components(g):
		j += 1
		if len(i) > max_len:
			max_len = len(i)
		
		if len(i) != 3172:
			sum_nodes = sum_nodes + len(i)
			
		nlist = list(i)
		if ('P69905' in nlist) or ('P68871' in nlist) :  
			sg = g.subgraph(nlist)
	
	# To plot the subgraph containing subAlpha and subBeta
	#nx.draw_networkx_nodes(sg, [P69905,P68871])
	#plt.show()
	
# CONSIDERATIONS ON GENERAL NETWORK
	
	# number of nodes and edges
	print("\nThe general network of direct interactions from IntAct, includes:")
	print('#NODES: ', g.number_of_nodes())
	print('#EDGES: ', g.number_of_edges())
	
	dgs = g.degree()
	ldgs = [(v,k) for k, v in dgs.items()]
	ldgs.sort()
	ldgs.reverse()
	print("\nProtein with the highest degree (",ldgs[0][0],") is ", ldgs[0][1])
	
	# average values of degrees, clustering and betweenness
	print('\nAverage value of degree: ', np.mean(list(dgs.values())))
	clst = nx.clustering(g)
	print('\nAverage clustering coefficient: ', np.mean(list(clst.values())))
	bwns = nx.betweenness_centrality_source(g, normalized = False)
	print('\nAverage betweenness : ', np.mean(list(bwns.values())))
	
	print("The network has ", j, "components")
	print("The biggest one has ", max_len, "nodes, and is the one containing alpha and beta subunits")
	print("Excluding the biggest subgraph, we have an average of ", sum_nodes/(j-1)," nodes")
	
	dgs_vals = list(dgs.values())
	hist,bins = np.histogram(dgs_vals, 50)
	#print (hist)
	#print(bins)
	
	#plot the histogram
	#plt.hist(dgs_vals, 50)
	#plt.title("A")
	#plt.xlabel("Degrees")
	#plt.ylabel("Number of Nodes")
	#plt.show()
	
	# Fit the distribution
	logx = []
	logy = []
	for i in range(len(hist)):
		if hist[i] == 0: continue
		logy.append(np.log10(hist[i]))
		logx.append(np.log10((bins[i+1]+bins[i])/2))
	#print (logx)
	#print (logy)
	reg = linregress(logx, logy)
	#print ("\nLINEAR REGRESSION")
	#print(reg)
	
	#create the fitted line
	yf = []
	for i in range(len(logx)):
		yf.append(logx[i]*reg[0]+reg[1])
	#plot linear regressiom
	#plt.plot(logx, logy, 'o')
	#plt.plot(logx, yf)
	#plt.title("B")
	#plt.xlabel("log(Mean Endpoints Degree Range)")
	#plt.ylabel("log(Number of Nodes)")
	#plt.show()

# CONSIDERATIONS ON COMPONENTS CONTAINING ALPHA AND BETA
# ALPHA = P69905, BETA = P68871	

	#DEGREES
	print ('sub ALPHA degree = ', sg.degree('P69905'))
	print( 'sub BETA degree = ', sg.degree('P68871'))
	
	#EDGES
	print ('sub ALPHA edges = ', sg.edges('P69905'))
	print( 'sub BETA edges = ', sg.edges('P68871'))
	
	# To plot a graph containing subAlpha and subBeta edges only
	#g_ab1 = nx.Graph()
	#g_ab2 = nx.Graph()
	#g_ab1.add_edges_from(sg.edges('P69905')+sg.edges('P68871'))
	#g_ab2.add_edges_from(sg.edges('P69905')+sg.edges('P68871')+sg.edges('P02008')+sg.edges('P0DMV8')+sg.edges('P02100')+sg.edges('P00387')+sg.edges('Q9NZD4')+sg.edges('P29474')+sg.edges('P69892'))
	#nx.draw(g_ab1, with_labels = True)
	#nx.draw(g_ab2, with_labels = True)	
	#plt.show()
	
	#CLUSTERING
	clust = nx.clustering(sg)
	print('sub ALPHA clustering = ', clust['P69905'])
	print('sub BETA clustering = ', clust['P68871'])
	
	#BETWEENNESS
	betweenness = nx.betweenness_centrality(sg, normalized = False)
	print('sub ALPHA betweenness = ', round(betweenness['P69905'], 3))
	print('sub BETA betweenness = ', round(betweenness['P68871'], 3))
	
	
#NETWORK ROBUSTNESS TEST
	# We remove node of subALPHA and check the change in betweenness
	# We do the same removing subBETA 
	
	nodes = list(sg.nodes())

	#remove subBETA
	nodes.remove('P68871')
	g_no_beta = sg.subgraph(nodes)
	betweenness1 = nx.betweenness_centrality(g_no_beta, normalized = False)

	nodes = list(sg.nodes())
	#remove subuALPHA
	nodes.remove('P69905')
	g_no_alpha = sg.subgraph(nodes)
	betweenness2 = nx.betweenness_centrality(g_no_alpha, normalized = False)
	
	print('Removing subunit beta, now the betweenness of sub alpha (P69905) is:')
	print(round(betweenness1['P69905'],3))
	print('Removing subunit alpha, now the betweenness of sub beta (P68871) is:')
	print(round(betweenness2['P68871'],3))


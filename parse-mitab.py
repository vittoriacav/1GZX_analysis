#!/usr/bin/env python3
import sys
import networkx as nx
import copy
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress

# "id1" and "id2" are the positions we are looking at. We want to look at
#  position 0 and 1 which are the swiss code (?)
def parse_mitab(filename, org, id1 = 0, id2 = 1):
	# define empty dictionary
	gi = {}
	# we open the file and read line by line as always
	f = open (filename, 'r')
	for line in f:
		if line[0] == '#': # that symbol indicate the header so we 
						   #  want to skip it
			continue
		line = line.rstrip()
		v = line.split('\t') # our file is tab separated
		# "v" is now a list where each item is a field of the line that 
		#  was tab-separated
		
		# Select only the lines in which both the interacting proteins 
		#  have a uniprotkb identifier
		# v[id1] = v[0] = unique identifier of interactor A (first column)
		# v[id2] = v[1] = unique identifier of interactor B (second column)
		
		# remember: ".find()" method finds th efirst occurrence of the 
		#  specified value. If the value is not found it returns -1
		
		if v[id1].find('uniprotkb:')==-1 or v[id2].find('uniprotkb:')==-1: 
			continue
		
		# We want consider only proteins interacting in a specific organism
		#  that in our case is human (taxid: 9606), hence we check the taxid
		#  for the 2 current interacting proteins
		# v[9] = NCBI taxonomy identifier for interactor A
		# v[10] = NCBI taxonomy identifier for interactor B
		
		if v[9].find('taxid:'+org)==-1 or v[10].find('taxid:'+org)==-1:
			continue
		
		# We split many times to obtain protein A and protein B names
		#  To better understand the process, check notes
		
		# In the end we will obtain for example:
		#  gid1 = 'P49418' (protein A)
		#  gid2 = 'O43426' (protein B)
		gid1 = v[id1].split('|')[0].split(':')[1].split('-')[0]
		gid2 = v[id2].split('|')[0].split(':')[1].split('-')[0]
		gids = [gid1, gid2]
		# following 2 steps onto identifiers are made
		#  to generate undirect edges later
		#  "list(gi.keys())" in the end (it's out return) will have pairs
		#  of unique interacting proteing, where each pair is sorted in 
		#  alphabetic order. This avoids duplicate interactions and also 
		#  duplicates in this sense: A-B = B-A
		#  To better understand the process, check notes
		gids.sort()
		gi[tuple(gids)] = True
		
		# example at iteration 8
		# list(gi.keys()) = [('O43426', 'Q99961'), ('P49418', 'Q05193'), ('O43426', 'P49418')]

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
	
#NOW WE CREATE THE GRAPH of subAlpha and subBeta

	g = nx.Graph()
	g.add_edges_from(gis)
	
	# To plot the whole graph based on gis
	#nx.draw(g)
	#plt.show()
	
	# We are using "nx.connected_components(g)" because it contains all
	#  components of our indirect graph. A graph in general can be composed
	#  of more than one components that are iduced subgraphs
	# We iterate all these "subgraphs" and check is subAlpha or subBeta 
	#  are present in the current subgraph; if it is present we create a
	#  variable with that subgraph which will be the subgraph that we will
	#  use for further analysis
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
			#print(len(nlist))  # --> 3172 is the number of connected nodes
								#     in the subgraph where subAlpha and 
								#     subBeta are present 
			sg = g.subgraph(nlist)
	
	# To plot the subgraph containing subAlpha and subBeta
	#nx.draw_networkx_nodes(sg, [P69905,P68871])
	#plt.show()
	
# GENERAL CONSIDERATIONS ON OUR NETWORK
	'''
	# number of nodes and edges
	print("\nThe general network of direct interactions from IntAct, includes:")
	print('#NODES: ', g.number_of_nodes())
	print('#EDGES: ', g.number_of_edges())
	
	# nodes (protein) with highest degree = protein which interact with 
	#  highest number of other proteins'''
	dgs = g.degree()
	'''ldgs = [(v,k) for k, v in dgs.items()]
	ldgs.sort()
	ldgs.reverse()
	print("\nProtein with the highest degree (",ldgs[0][0],") is ", ldgs[0][1])
	
	# average values of degrees, clustering and betweenness
	print('\nAverage value of degree: ', np.mean(list(dgs.values())))
	clst = nx.clustering(g)
	print('\nAverage clustering coefficient: ', np.mean(list(clst.values())))
	bwns = nx.betweenness_centrality_source(g, normalized = False)
	print('\nAverage betweenness : ', np.mean(list(bwns.values())))
	
	print("The network has ", j, "components (subgraphs)")
	print("The biggest one has ", max_len, "nodes, and is the one containing alpha and beta subunits")
	print("Excluding the biggest subgraph, we have an average of ", sum_nodes/(j-1)," nodes")
	
	# retrieve the list of values for the degree of our nodes'''
	dgs_vals = list(dgs.values())
	plt.hist(dgs_vals, 50)
	plt.title("A")
	plt.xlabel("Degrees")
	plt.ylabel("Number of Nodes")
	plt.show()
	'''
	
	# numpy function that create the histogram and return it as a list 
	#  1st argument: list of values
	#  2nd argument: number of bins you want to divide the data
	# return : a tuple --> first arg is the list of number of occurrences  
	#					   (#nodes)in each bin
    #				   --> second arg is the list of range of the bins 
	
	hist,bins = np.histogram(dgs_vals, 50)
	#print (hist)
	#print(bins)

	# We want now to fit this distribution
	# We want to use linear regression but we have a power law function
	#  so we reduce our function to a linear one using the logarithm
	#  y = b*x^(-a)   -------------->  log(y) = -a*log(x) + log(b)
	# So we need to transform our list of numbers in list of logarithms
	# In our cast "hist" are our y and "bins" are our x
	
	logx = []
	logy = []
	for i in range(len(hist)):
		# Be careful because log(0) is not possible
		if hist[i] == 0: continue
		logy.append(np.log10(hist[i]))
		# Note that the range of x is equal to the range of y but +1
		#  because for each value of the hist we have starting and ending 
		#  point (and for next value of hist the starting point is the ending 
		#  point of the previous hist value)
		#  So we can take and average between starting and ending point
		logx.append(np.log10((bins[i+1]+bins[i])/2))
	#print (logx)
	#print (logy)
	reg = linregress(logx, logy)
	#print ("\nLINEAR REGRESSION")
	#print(reg)
	
	#reg[0] = slope,  reg[1] = intercept ---> we can draw the fitted line
	yf = []
	for i in range(len(logx)):
		yf.append(logx[i]*reg[0]+reg[1])
	print(reg)
	plt.plot(logx, logy, 'o')
	plt.plot(logx, yf)
	plt.title("B")
	plt.xlabel("log(Mean Endpoints Degree Range)")
	plt.ylabel("log(Number of Nodes)")
	plt.show()
	'''
	'''	
# ALPHA = P69905, BETA = P68871	

	#TASK 4: network topology
	
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
	
	
	# TASK 5: network robustness
	# We remove node of subALPHA and check the change in betweenness
	# We do the same remoivng subBETA this time
	
	# REMEMBER BETWEENNESS CENTRALITY: is a measure of centrality in a 
	#  graph based on shortest path. For every pair of vertices in a 
	#  connected graph, there exists at least one shortest path between
	#  them. The betweenness centrality fot each vertes is the number of
	#  of these shortest paths that pass through the vertex
	
	nodes = list(sg.nodes())
	
	#remove subBETA
	nodes.remove('P68871')
	g_no_beta = sg.subgraph(nodes)
	#deepcopy alternative for previous commands
	#g_no_beta = copy.deepcopy(sg)
	#g_no_beta.remove_node('P68871')
	betweenness1 = nx.betweenness_centrality(g_no_beta, normalized = False)

	nodes = list(sg.nodes())
	#remove subuALPHA
	nodes.remove('P69905')
	g_no_alpha = sg.subgraph(nodes)
	#deepcopy alternative for previous commands
	#g_no_alpha = copy.deepcopy(sg)
	#g_no_alpha.remove_node('P69905')
	betweenness2 = nx.betweenness_centrality(g_no_alpha, normalized = False)
	
	print('Removing subunit beta, now the betweenness of sub alpha (P69905) is:')
	print(round(betweenness1['P69905'],3))
	print('Removing subunit alpha, now the betweenness of sub beta (P68871) is:')
	print(round(betweenness2['P68871'],3))
	
	#As we expected, removing subunit alpha, that is the main node connecting 
	# to beta, the betweenness centrality of beta goes to 0
	'''

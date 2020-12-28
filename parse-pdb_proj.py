#!/usr/bin/env python3
import sys
import numpy as np

# Function parsing the PDB file
def parse_pdb(filename, chain):
	pdb_coords={}
	f=open(filename,'r')
	for line in f:
		line=line.rstrip()
		# We would like to retain lines containing:
		#  - the atoms of the residues of each chain of 
		#    hemoglobin --> lines that start with ATOM
		#  - the atoms of the residues of both heme group (HEM) 
		#    and oxygen group (OXY) --> lines that start with HETATM
		
		# Check the field ATOM and HETATM at the beginning of the line
		# Check the correct chain in column 21
		
		# As we can see from the 2 examples of possible lines from the
		#  pdb file (one the ATOM and one for the HETATM), the chain is
		#  always specified at column 21
		
		#                     21    
		#					  |
		#ATOM      1  N   VAL A   1      18.432  -2.931   3.579  1.00 37.68           N
		#HETATM 4389  CHA HEM A1142      17.735  18.786  16.239  1.00 27.82
		            
		if (line[:4] == 'ATOM' or line[:6] == 'HETATM') and line[21] == chain: 
			pass
		else:
			continue
			
		#NUMBER OF THE RESIDUE
		#                      22 25   
		#					   |  |    
		#ATOM      1  N   VAL A   1      18.432  -2.931   3.579  1.00 37.68           N
		#HETATM 4389  CHA HEM A1142      17.735  18.786  16.239  1.00 27.82	
		
		resn = int(line[22:26].strip()) #22:26 takes the contents of from col 22 to col 25, excluding col 26
		
		#NAME OF THE RESIDUE
		#                 1719                         
		#				  | |	         
		#ATOM      1  N   VAL A   1      18.432  -2.931   3.579  1.00 37.68           N
		#HETATM 4389  CHA HEM A1142      17.735  18.786  16.239  1.00 27.82	
		
		res = line[17:20].strip()
		
		# In HETATM lines we also find HOH (water), used to produce the crystal
		#  structure. We need to remove them since we are not interested in them
		if res == 'HOH': continue
		
		# We extract now the COORDINATES of the atoms for both ATOM and HETATM

		atom = line[12:16].strip()
		x = float(line[30:38])
		y = float(line[38:46])
		z = float(line[46:54])
		coord = [x,y,z]
		
		# Initialize the dictionary with the number of each residue: resn
		# and "attach" to each key an empty dictionary, that we will use later
		pdb_coords[resn] = pdb_coords.get(resn,{}) #if this key already exists bc we created it with previous line
												   # is will not add another one but will use the exixsting one for next step
		# Add the atom's coordinates and the name of the residue we are considering
		pdb_coords[resn][atom] = coord
		pdb_coords[resn]['name'] = res
		
		#pdb_coord DICTIONARY STRUCTURE
		
		# example
		# {1: {'CB': [20.659, -3.754, 2.825], 'name': 'VAL', 
		# 'CA': [19.662, -2.549, 2.806], 'CG2': [21.982, -3.272, 2.245], 
		# 'CG1': [20.109, -4.992, 2.222], 'N': [18.432, -2.931, 3.579], 
		# 'O': [18.421, -2.497, 0.695], 'C': [19.282, -1.939, 1.441]}}

		# We have a dictionary where keys are the residue numbers. In the
		#  example we are considering residue number 1 which is a valine 
		#  (VAL) and in the pdb file is represented by these lines:
		
		#ATOM      1  N   VAL A   1      18.432  -2.931   3.579  1.00 37.68           N  
		#ATOM      2  CA  VAL A   1      19.662  -2.549   2.806  1.00 35.41           C  
		#ATOM      3  C   VAL A   1      19.282  -1.939   1.441  1.00 34.04           C  
		#ATOM      4  O   VAL A   1      18.421  -2.497   0.695  1.00 33.95           O  
		#ATOM      5  CB  VAL A   1      20.659  -3.754   2.825  1.00 35.59           C  
		#ATOM      6  CG1 VAL A   1      20.109  -4.992   2.222  1.00 37.84           C  
		#ATOM      7  CG2 VAL A   1      21.982  -3.272   2.245  1.00 36.73           C  
		
		# The argument of the key is another dictionary where keys are 
		#  the names of the atoms composing the residue (in this case
		#  atoms composing valine, residue 1).
		#  An additional special key is also present, called 'name', that
		#  indicates the name of the residue we are considering (in this
		#  case, VAL)
		
	return pdb_coords


# Function calculating the distance between two points
def get_distance(coord1,coord2):
	return np.sqrt((coord1[0]-coord2[0])**2+\
				   (coord1[1]-coord2[1])**2+\
                   (coord1[2]-coord2[2])**2)


# Function returning the distances below a certain threshold
# we want to calculate the distance between 2 residues
def get_min_dist(res_coords1, res_coords2, th = 3.5):
	atm_dists = []
	# these lists contain all the atom names of each one of the 2 residues
	keys1 = list(res_coords1.keys())
	keys1.remove('name')
	keys2 = list(res_coords2.keys())
	keys2.remove('name')
	
	for k1 in keys1:
		for k2 in keys2:
			dist = get_distance(res_coords1[k1], res_coords2[k2])
			if dist <= th: 
				atm_dists.append([res_coords1['name'], k1, res_coords2['name'], k2, dist])
	# we want to return a list of atoms that have a distance below a fixed threshold
	return atm_dists

def get_het_dist(pdb_coords, heme, oxy):
	# We create a list called "keys" that will contain all the number of 
	#  the residues of chain A/B/C/D. We remove infact residue numbers 
	#  corresponding to the heme and oxygen groups
	keys = list(pdb_coords.keys())
	keys.remove(heme)
	keys.remove(oxy)
	keys.sort()	
	# Header to make output more comprehensible
	print('res\tres\tatom\tgroup\tatom\tdist')
	print('numb\tname\tres\tname\taroup\tres-group\n')

	# For each residue number in the list of residue numbers...
	for k in keys:
		# ... use "get_min_dist" to retrieve the distance between 2
		#     residues (one from the chain A/B/C/D and the HEM/OXY
		#     group, which is one residue itself). 
		#     "pdb_coords[heme/oxy]" will always be the same 
		
		# example of possible return of "get_min_dist" function: the
		#  in order:
		#	- name of the residue belgoning to chain A/B/C/D
		#	- specific atom of that residue
		#	- HEME or OXY depending on the residue/group we are considering
		#	- specific atom of the heme/oxy group/residue
		#	- distance between the 2 atoms 
		
		# NOTE: only distances below 3.5 A (threshold) will be returned
		
		# example: in this case (k = 42) atom O of TYR residue (42) is
		#  distance 3.341 from atom CMD of HEM group/residue
		# [['TYR', 'O', 'HEM', 'CMD', 3.3414040462057275]]
		atm_dists = get_min_dist(pdb_coords[k], pdb_coords[heme])
		if len(atm_dists) > 0:
			'''print("############################################à")
			print("K = ", k)
			print("HEME")'''
			for l in atm_dists:
				print(str(k)+'\t'+'\t'.join([str(j) for j in l]))
		atm_dists = get_min_dist(pdb_coords[k], pdb_coords[oxy])
		if len(atm_dists) > 0:
			'''print("OXY")'''
			for l in atm_dists:
				print(str(k)+'\t'+'\t'.join([str(j) for j in l]))
				
def get_chain_dist(pdb_coords1, pdb_coords2, chain1, chain2):
	keys1 = list(pdb_coords1.keys())
	keys2 = list(pdb_coords2.keys())
	keys1.sort()
	keys2.sort()
	# Header to make output more comprehensible
	print("\nnumb\tnumb\tres1\tatom\tres2\tatom\tdist")
	print("res1\tres2\tname\tres1\tname\tres2\tre1-res2\n")
	for k1 in keys1:
		name1 = pdb_coords1[k1]['name']
		if name1 == 'OXY' or name1 == 'HEM': continue
		for k2 in keys2:
			name2 = pdb_coords2[k2]['name']
			if name2 == 'OXY' or name2 == 'HEM': continue
			atm_dists = get_min_dist(pdb_coords1[k1], pdb_coords2[k2], 4.0)
			if len(atm_dists) > 0:
				for l in atm_dists:
					print(chain1+str(k1)+'\t'+chain2+str(k2)+'\t'+'\t'.join([str(j) for j in l]))
			# OUTPUT STRUCTURE
			# In order:
			#	- chain 1 name + considered residue number belonging to chain 1
			#	- chain 2 name + onsidered residue number belonging to chain 2
			#	- name of the residue belgoning to chain 1
			#	- specific atom of residue belgoning to chain 1
			#	- name of the residue belgoning to chain 2
			#	- specific atom of residue belgoning to chain 2
			#	- distance between the 2 atoms 
			
			
			
			
if __name__ == '__main__':
	# As input we have the file name of the pdb structure, in our case: 1gzx.pdb
	filename=sys.argv[1]
	
	# We obtain a dictionary for each chain where all residues are present
	#  along with the coordinates of the atoms composing each residue. Of 
	#  course also the heme and oxygen groups are present in the dictionary 
	chain= ['A', 'B', 'C', 'D']
	pdb_coordsA=parse_pdb(filename,chain[0])
	pdb_coordsB=parse_pdb(filename,chain[1])
	pdb_coordsC=parse_pdb(filename,chain[2])
	pdb_coordsD=parse_pdb(filename,chain[3])
	
	
	# TASK 1
	# Considering a minimum of 3.5 A, identify the residues, belonging to
	#  hemoglobin subunits, in proximity of the heme and oxygen group
	
	# HEME   OXY  CHAIN
	# resn   resn   |
	#  |      |     |
	# 1142, 1143    A
	# 1290, 1291    B
	# 1542, 1543    C
	# 1690, 1691    D
	
	'''print("\nInteracting residues between CHAIN A and HEM/OXY groups\n")
	#1142 is the residue number of the HEM group
	#1143 is the residue number of the OXY group
	get_het_dist(pdb_coordsA, 1142, 1143)
	print("\nInteracting residues between CHAIN B and HEM/OXY groups\n")
	get_het_dist(pdb_coordsB, 1290, 1291)
	print("\nInteracting residues between CHAIN C and HEM/OXY groups\n")
	get_het_dist(pdb_coordsC, 1542, 1543)
	print("\nInteracting residues between CHAIN D and HEM/OXY groups\n")
	get_het_dist(pdb_coordsD, 1690, 1691)'''
	
	
	# TASK 2
	# Using the same procedure of TASK 1, identify the possible 
	#  interacting residues between each pair of monomers
	
	'''
	#Chain A and Chain B
	print("\nInteracting residues between CHAIN ",chain[0], " and CHAIN ", chain[1],"\n")
	get_chain_dist(pdb_coordsB, pdb_coordsA, chain[1], chain[0])
	# Chain B and Chain C
	print("\nInteracting residues between CHAIN ",chain[1], " and CHAIN ", chain[2],"\n")
	get_chain_dist(pdb_coordsB, pdb_coordsC, chain[1], chain[2])
	# Chain C and Chain D
	print("\nInteracting residues between CHAIN ",chain[2], " and CHAIN ", chain[3],"\n")
	get_chain_dist(pdb_coordsC, pdb_coordsD, chain[2], chain[3])'''
	# Chain A and Chain D
	print("\nInteracting residues between CHAIN ",chain[3], " and CHAIN ", chain[0],"\n")
	get_chain_dist(pdb_coordsD, pdb_coordsA, chain[3], chain[0])
	'''# Chain A and Chain C
	print("\nInteracting residues between CHAIN ",chain[0], " and CHAIN ", chain[2],"\n")
	get_chain_dist(pdb_coordsA, pdb_coordsC, chain[0], chain[2])
	# Chain B and Chain D
	print("\nInteracting residues between CHAIN ",chain[1], " and CHAIN ", chain[3],"\n")
	get_chain_dist(pdb_coordsB, pdb_coordsD, chain[1], chain[3])'''
  

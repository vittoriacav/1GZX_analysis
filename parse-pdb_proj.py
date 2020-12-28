#!/usr/bin/env python3
import sys
import numpy as np

# Function parsing the PDB file
def parse_pdb(filename, chain):
	pdb_coords={}
	f=open(filename,'r')
	for line in f:
		line=line.rstrip()		            
		if (line[:4] == 'ATOM' or line[:6] == 'HETATM') and line[21] == chain: 
			pass
		else:
			continue
			
		resn = int(line[22:26].strip()) 
		res = line[17:20].strip()
		
		# In HETATM lines we also find HOH (water), used to produce the crystal
		#  structure. We need to remove them
		if res == 'HOH': continue
		
		# We extract now the COORDINATES of the atoms for both ATOM and HETATM

		atom = line[12:16].strip()
		x = float(line[30:38])
		y = float(line[38:46])
		z = float(line[46:54])
		coord = [x,y,z]
		
		pdb_coords[resn] = pdb_coords.get(resn,{}) 
		pdb_coords[resn][atom] = coord
		pdb_coords[resn]['name'] = res		
	return pdb_coords


# Function calculating the distance between two points
def get_distance(coord1,coord2):
	return np.sqrt((coord1[0]-coord2[0])**2+\
				   (coord1[1]-coord2[1])**2+\
                   (coord1[2]-coord2[2])**2)


# Function returning the distances below a certain threshold
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
	return atm_dists

def get_het_dist(pdb_coords, heme, oxy):
	keys = list(pdb_coords.keys())
	keys.remove(heme)
	keys.remove(oxy)
	keys.sort()	
	# Header to make output more comprehensible
	print('res\tres\tatom\tgroup\tatom\tdist')
	print('numb\tname\tres\tname\taroup\tres-group\n')
	
	for k in keys:
		
		atm_dists = get_min_dist(pdb_coords[k], pdb_coords[heme])
		if len(atm_dists) > 0:
			for l in atm_dists:
				print(str(k)+'\t'+'\t'.join([str(j) for j in l]))
		
		atm_dists = get_min_dist(pdb_coords[k], pdb_coords[oxy])
		if len(atm_dists) > 0:
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
			
if __name__ == '__main__':
	filename=sys.argv[1]

	chain= ['A', 'B', 'C', 'D']
	pdb_coordsA=parse_pdb(filename,chain[0])
	pdb_coordsB=parse_pdb(filename,chain[1])
	pdb_coordsC=parse_pdb(filename,chain[2])
	pdb_coordsD=parse_pdb(filename,chain[3])
	
	
	# Identify the residues, belonging to hemoglobin subunits, in proximity of the heme and oxygen group
	print("\nInteracting residues between CHAIN A and HEM/OXY groups\n")
	get_het_dist(pdb_coordsA, 1142, 1143)
	print("\nInteracting residues between CHAIN B and HEM/OXY groups\n")
	get_het_dist(pdb_coordsB, 1290, 1291)
	print("\nInteracting residues between CHAIN C and HEM/OXY groups\n")
	get_het_dist(pdb_coordsC, 1542, 1543)
	print("\nInteracting residues between CHAIN D and HEM/OXY groups\n")
	get_het_dist(pdb_coordsD, 1690, 1691)
	
	
	# Identify the possible interacting residues between each pair of monomers
	
	
	#Chain A and Chain B
	print("\nInteracting residues between CHAIN ",chain[0], " and CHAIN ", chain[1],"\n")
	get_chain_dist(pdb_coordsB, pdb_coordsA, chain[1], chain[0])
	# Chain B and Chain C
	print("\nInteracting residues between CHAIN ",chain[1], " and CHAIN ", chain[2],"\n")
	get_chain_dist(pdb_coordsB, pdb_coordsC, chain[1], chain[2])
	# Chain C and Chain D
	print("\nInteracting residues between CHAIN ",chain[2], " and CHAIN ", chain[3],"\n")
	get_chain_dist(pdb_coordsC, pdb_coordsD, chain[2], chain[3])
	# Chain A and Chain D
	print("\nInteracting residues between CHAIN ",chain[3], " and CHAIN ", chain[0],"\n")
	get_chain_dist(pdb_coordsD, pdb_coordsA, chain[3], chain[0])
	# Chain A and Chain C
	print("\nInteracting residues between CHAIN ",chain[0], " and CHAIN ", chain[2],"\n")
	get_chain_dist(pdb_coordsA, pdb_coordsC, chain[0], chain[2])
	# Chain B and Chain D
	print("\nInteracting residues between CHAIN ",chain[1], " and CHAIN ", chain[3],"\n")
	get_chain_dist(pdb_coordsB, pdb_coordsD, chain[1], chain[3])
  

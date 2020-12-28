#!/usr/bin/env python3
import sys
import numpy as np

#solvent AA accessibility = saa
saa = {'A': 115, 'L':170, 'R':225, 'K':200, 'N':160, 'M':185, 'D':150, 'F':210,\
 'C': 135, 'P': 145, 'Q': 180, 'S': 115, 'E':190, 'T':140, 'G':75, 'W':255, 'H':195,
 'Y':230, 'I':175, 'V':155}
		
# Function parsing the DSSP file
def parse_dssp(filename,chain):
	dssp_list=[]
	c=0
	f=open(filename)
	for line in f:
		line=line.rstrip()
		if line.find('  #  RESIDUE')>-1: 
			c=1 
			continue
		# Check if the state variable is 1
		# Select the correct chain
		if c==0 or line[11]!=chain: continue
		
		aa=line[13]
		if aa.islower(): aa='C'
		if aa=='!': continue
		resn=int(line[5:10])
		ss=line[16]
		if ss==' ': ss='C'
		acc=float(line[34:38])
		phi=float(line[103:109])
		psi=float(line[109:115])
		if saa.get(aa,0) == 0: continue
		rsa = round(acc/saa[aa], 3)
		rsa = min(rsa,1.0)
		aa_dssp=[resn, aa, ss, acc, rsa]
		dssp_list.append(aa_dssp)
	return dssp_list

def get_chain_acc(dssp_list):
	accs=[]
	for aa in dssp_list:
		accs.append(aa[3])
	return accs

if __name__ == '__main__':
	filename=sys.argv[1]
	chain=sys.argv[2]
	dssp_list=parse_dssp(filename,chain)
	#Uncomment to get normlized solvent accessibility
	#print ('Total normalized solvent accessiblity for chain ',chain,':', sum(accs))
	
	accs = get_chain_acc(dssp_list)
	for i_dssp in dssp_list:
		print ('\t'.join([str(i) for i in i_dssp]))
	

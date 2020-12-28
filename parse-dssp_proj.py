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
	# Initialize state variable c=0
	c=0
	f=open(filename)
	for line in f:
		line=line.rstrip()
		# Check for the beginning of data
		# Change the state variable c=1
		if line.find('  #  RESIDUE')>-1: 
			c=1 #finding that word we don't want to read that line
			continue
		# Check if the state variable is 1
		# Select the correct chain
		if c==0 or line[11]!=chain: continue
		# Read all the variables resn, aa, acc, phi, psi
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
		
'''
# Get sequence and secondary structure
def get_ss(dssp_list):
	seq=''
	ss=''
	for aa in dssp_list:
		ss=ss+aa[2]
		seq=seq+aa[1]
	return seq,ss
	

# Count residues with a given secondary structure
def count_ss(ss,ss_type):
	c=0
	for i in ss:
		if i==ss_type: c=c+1
	return c 


# Return the list of all phi and psi angles
def get_ss_angle(dssp_list,ss_type):
	angles=[]
	for aa in dssp_list:
		if aa[2]==ss_type: angles.append((aa[-2],aa[-1]))
	return angles

# Return the list of residues accessibility	
def get_aa_acc(dssp_list,aa_type):
	accs=[]
	for aa in dssp_list:
		if aa[1]==aa_type: accs.append(aa[3])
	return accs
'''

def get_chain_acc(dssp_list):
	accs=[]
	for aa in dssp_list:
		accs.append(aa[3])
	return accs

if __name__ == '__main__':
	filename=sys.argv[1]
	chain=sys.argv[2]
	dssp_list=parse_dssp(filename,chain)
	# "dssp_list" is a list of lists
	# Each list contains:
	#	- residue number
	#	- residue one-letter-name
	#	- secondary structure type:
	#			* H = Alpha helix
	#			* B = Isolated beta-bridge residue
	#			* E = Strand
	#			* G = 3-10 helix
	#			* I = Pi helix
	#			* T = Turn
	#			* S = Bend
	#			* - = None
	#	- solvent accessibility: number of water molecules in contact with
	#	  this residue *10. or residue water exposed surface in Armstrong**2 
	#	- normalized solvent accessibility
	
	#example of one of the list inside the main list:
	#[1, 'V', 'C', 157.0, 1.0]
	
	accs = get_chain_acc(dssp_list)
	for i_dssp in dssp_list:
		print ('\t'.join([str(i) for i in i_dssp]))
	#print ('Total normalized solvent accessiblity for chain ',chain,':', sum(accs))

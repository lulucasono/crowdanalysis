# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 10:21:49 2018

@author: wth
"""

import pandas as pd
#import numpy as np
import sys
import ast

input_file_name = sys.argv[1]
output_file_name = sys.argv[2]


input_seq = pd.read_csv(input_file_name,header=None, encoding = "UTF-8", sep='\t')

max_id_of_people = max(input_seq[0])
tth = 1800

class Node:
	def __init__(self, cluster, i, j):
		self.cluster = cluster
		self.ipos = i
		self.jpos = j
	
	def __str__(self):
		return str(self.cluster)+"["+str(self.ipos)+","+str(self.jpos)+"]"
	
	def __eq__(self, other):
		"""Overrides the default implementation"""
		if isinstance(other, Node):
			return (self.cluster == other.cluster) & (self.ipos == other.ipos) & (self.jpos == other.jpos)
		return False
	
	def __repr__(self):
		return str(self)
	
	def __hash__(self):
		return hash(str(self))

	
def ExtendSequence(users_seq,seq,i,j,tth):
	res = []
	seqi = ast.literal_eval(users_seq.iloc[i,1])
	seqj = ast.literal_eval(users_seq.iloc[j,1])
	timei = ast.literal_eval(users_seq.iloc[i,3])
	timej = ast.literal_eval(users_seq.iloc[j,3])
	lastNode = seq[-1]
	for ii in range(lastNode.ipos+1,len(seqi)):
		deltai = sum(timei[lastNode.ipos:ii])
		for jj in range(lastNode.jpos+1,len(seqj)):
			deltaj = sum(timej[lastNode.jpos:jj])
			if (seqi[ii]==seqj[jj]) & (abs(deltai-deltaj)<=tth):
				res.append(seq.copy())
				res[-1].append(Node(seqi[ii],ii,jj))
				
	return res

class LCS:
	def __init__(self, timei, timej, tth):
		self.timei = timei
		self.timej = timej
		self.tth = tth
		self.lcss = []
		
	def getAllLcs(self, flag, seqi, lcs, maxSublen, i, j):
		while (i>0)&(j>0)&(len(lcs)<maxSublen):
			direction = flag[i][j]
			if direction == 2:
				lcs.append(Node(seqi[i-1],i-1,j-1))
				i -= 1
				j -= 1
			else:
				if direction == 1:
					i -= 1
				elif direction == 3:
					j -= 1
				else:
					self.getAllLcs(flag, seqi, lcs, maxSublen, i-1, j)
					self.getAllLcs(flag, seqi, lcs, maxSublen, i, j-1)
		if len(lcs)==maxSublen:
			self.lcss.append(lcs[::-1])
		

def SequenceMatching(users_seq, i, j,tth):
	seqi = ast.literal_eval(users_seq.iloc[i,1])
#	print(seqi)
	seqj = ast.literal_eval(users_seq.iloc[j,1])
#	print(seqj)
	timei = ast.literal_eval(users_seq.iloc[i,3])
	timej = ast.literal_eval(users_seq.iloc[j,3])
	
	leni=len(seqi)
	lenj=len(seqj)
	c=[[0 for x in range(lenj+1)] for x in range(leni+1)]
	flag=[[0 for x in range(lenj+1)] for x in range(leni+1)]
	for x in range(leni):
		for y in range(lenj):
			if seqi[x]==seqj[y]:
				c[x+1][y+1]=c[x][y]+1
				flag[x+1][y+1] = 2
			elif c[x][y+1]>c[x+1][y]:
				c[x+1][y+1]=c[x][y+1]
				flag[x+1][y+1] = 1
			elif c[x][y+1]<c[x+1][y]:
				c[x+1][y+1]=c[x+1][y]
				flag[x+1][y+1] = 3
			else:
				c[x+1][y+1]=c[x][y+1]
				flag[x+1][y+1] = 4
				
	maxStep = max(max(c))
	
	ob_LCS = LCS(timei, timej, tth)
	lcs = []
	lcss = []
	for x in range(leni,0,-1):
		for y in range(lenj,0,-1):
			if c[x][y] >= maxStep:
				lcs = []
				ob_LCS.getAllLcs(flag, seqi, lcs, maxStep, x, y)
	lcs_set = set(tuple(x) for x in ob_LCS.lcss)
	lcss = [ list(x) for x in lcs_set ]
	
	return lcss

def SimilarityMeasureOfSequence(seq,i,j):
	seqNi = ast.literal_eval(input_seq.iloc[i,2])
	seqNj = ast.literal_eval(input_seq.iloc[j,2])
	timei = ast.literal_eval(input_seq.iloc[i,3])
	timej = ast.literal_eval(input_seq.iloc[j,3])
	if len(seq)==0:
		return 0
	sm = []
	res = 0
#	n = 0
	x = 1
	lastx = x-1
#	while x<len(seq):
#		node = seq[x]
#		n += 1
#		if x>0:
#			deltai = sum(timei[seq[x-1].ipos:seq[x].ipos])
#			deltaj = sum(timej[seq[x-1].jpos:seq[x].jpos])
#			if abs(deltai-deltaj)<=tth:
#				sm += 2**(n-1)*min(seqNi[node.ipos],seqNj[node.jpos])
#			else:
#				n = 1
#				sm += 2**(n-1)*min(seqNi[node.ipos],seqNj[node.jpos])
#		else:
#			sm += min(seqNi[node.ipos],seqNj[node.jpos])
#	sm = 2**(m-1) * sm	
#	for x in range(1,len(seq)):
#		deltai = sum(timei[seq[x-1].ipos:seq[x].ipos])
#		deltaj = sum(timej[seq[x-1].jpos:seq[x].jpos])
#		if abs(deltai-deltaj)<=tth:
#			sm *= 2
	for node in seq:
		sm.append( min(seqNi[node.ipos],seqNj[node.jpos]))
	while x < len(seq):
		deltai = sum(timei[seq[x-1].ipos:seq[x].ipos])
		deltaj = sum(timej[seq[x-1].jpos:seq[x].jpos])
		if abs(deltai-deltaj)>tth:
			res += 2**(x-lastx)*sum(sm[lastx:x+1])
			lastx = x-1
		x += 1
		
	if lastx == 0:
		res += 2**(x+1-lastx)*sum(sm[lastx:x+1])
	
	return res

def SimilarityOfLayer(seq_list,i,j):
#	print("Enter SimilarityOfLayer")
	ni = input_seq.iloc[i,4]
	nj = input_seq.iloc[j,4]
	s = 0;
	for seq in seq_list:
		s += SimilarityMeasureOfSequence(seq,i,j)
	s = s/(ni*nj)
	return s

similarity = [[0 for i in range(max_id_of_people+1)] for i in range(max_id_of_people+1)]

for j in range(len(input_seq[0])):
	print("=============================================")
	print("Number of person:",j)
	pidj = input_seq.iloc[j,0]
	print("perseon id:",	pidj)
	for i in range(j+1,len(input_seq[0])):
		pid = input_seq.iloc[i,0]
		if not((pidj==4)&(pid==59)|(pid==4)&(pidj==59)):
			max_len_seq_table = SequenceMatching(input_seq,j,i,tth)
#			print(max_len_seq_table)
			similarity[pidj][pid] = SimilarityOfLayer(max_len_seq_table, j, i)
			similarity[pid][pidj] = SimilarityOfLayer(max_len_seq_table, j, i)
			print("Similarity("+str(pidj)+","+str(pid)+"):",similarity[pidj][pid],(similarity[pid][pidj]==len(max_len_seq_table))|(len(max_len_seq_table)>0))
#			print()


df = pd.DataFrame(similarity)
df.to_csv(output_file_name,header=None, index=None, sep='\t')

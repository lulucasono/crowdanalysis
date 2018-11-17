# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 08:20:24 2018

@author: wth
"""

import pandas as pd
#import numpy as np
import sys
import ast

input_file_name = sys.argv[1]
output_file_name = sys.argv[2]
mlength = sys.argv[3]

input_seq = pd.read_csv(input_file_name,header=None, encoding = "UTF-8", sep='\t')
#input_seq = pd.read_csv("users_Seq_001",header=None, encoding = "UTF-8", sep='\t')
#input_seq2 = pd.read_csv("users_Seq_01",header=None, encoding = "UTF-8", sep='\t')
#input_seq3 = pd.read_csv("users_Seq_1",header=None, encoding = "UTF-8", sep='\t')

#person_idx = input_seq.index[input_seq[0]==int(person_id)][0]
#print(idx_001)
#idx_01 = input_seq2.index[input_seq2[0]==int(person_id)][0]
#print(idx_01)
#idx_1 = input_seq3.index[input_seq3[0]==int(person_id)][0]
#print(idx_1)

max_id_of_people = max(input_seq[0])
#max_id_of_people_01 = max(input_seq2[0])
#max_id_of_people_1 = max(input_seq3[0])

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
				if(len(res)==0):
					res = [seq.copy()]
				res[-1].append(Node(seqi[ii],ii,jj))
				res.append(seq.copy())
					
	return res

def SequenceMatching(users_seq, i,j,maxLength,tth):
	res = []
	step = 0
	seqi = ast.literal_eval(users_seq.iloc[i,1])
	seqj = ast.literal_eval(users_seq.iloc[j,1])

	for x in range(len(seqi)):
		for y in range(len(seqj)):
			if (seqi[x]==seqj[y]):
				node = Node(seqi[x],x,y)
				res.append([node])
	
	step=1
#	print(res)
	while step<maxLength:
		extended_seqs = []
		for s in res:
			if len(s)==step:
				extended_seqs.append(ExtendSequence(users_seq,s,i,j,tth))	

		
		for es in extended_seqs:
			res = res + es
		
		if len(res)>0:
			m = max(len(l) for l in res)
		else:
			break
		
		step += 1
		res = [e for e in res if len(e)==m]
		
	return res

def SimilarityMeasureOfSequence(seq,i,j):
	seqNi = ast.literal_eval(input_seq.iloc[i,2])
	seqNj = ast.literal_eval(input_seq.iloc[j,2])
	if len(seq)==0:
		return 0
	m = len(seq)
	sm = 0
	for node in seq:
		sm += min(seqNi[node.ipos],seqNj[node.jpos])
	sm = 2**(m-1) * sm	
	
	return sm

def SimilarityOfLayer(seq_list, i, j):
	ni = input_seq.iloc[i,4]
	nj = input_seq.iloc[j,4]
	s = 0;
	for seq in seq_list[input_seq.iloc[i,0]][input_seq.iloc[j,0]]:
		s += SimilarityMeasureOfSequence(seq,i,j)
	s = s/(ni*nj)
	return s

#seqi = ast.literal_eval(input_seq.iloc[1,1])
#seqj = ast.literal_eval(input_seq.iloc[1,1])
#
#res = SequenceMatching(1,1,5,10800)
#print("seqi:",seqi)
#print("seqj:",seqj)
#print("res:")
#for r in res:
#	print(r)

tth = 10800
max_length = int(mlength)
#max_len_seq_table = [[[] for i in range(max_id_of_people+1)] for i in range(max_id_of_people+1)]
#seqi = ast.literal_eval(input_seq.iloc[int(person_id),1])
#print("seq_"+str(input_seq.iloc[int(person_id),0])+":",seqi)
max_len_seq_table = [[-1 for i in range(max_id_of_people+1)] for i in range(max_id_of_people+1)]
similarity = [[0 for i in range(max_id_of_people+1)] for i in range(max_id_of_people+1)]
for j in range(len(input_seq[0])):
	print("=============================================")
	pidj = input_seq.iloc[j,0]
	print("perseon id:",	pidj)
	for i in range(j,len(input_seq[0])):
		pid = input_seq.iloc[i,0]
		max_len_seq_table[pidj][pid] = SequenceMatching(input_seq,j,i,max_length,tth)
		similarity[pidj][pid] = SimilarityOfLayer(max_len_seq_table, j, i)
		similarity[pid][pidj] = SimilarityOfLayer(max_len_seq_table, j, i)
#		print(pidj,"-",pid,":",max_len_seq_table[input_seq.iloc[i,0]])
		print("Similarity("+str(pidj)+","+str(pid)+"):",similarity[pidj][pid])

#for i in range(len(input_seq[0])):
#	print(i)
#	for j in range(i,len(input_seq[0])):
#		max_len_seq_table[input_seq.iloc[i,0]][input_seq.iloc[j,0]] = SequenceMatching(i,j,max_length,tth)
#		print(input_seq.iloc[i,0],"-",input_seq.iloc[j,0],":",max_len_seq_table[input_seq.iloc[i,0]][input_seq.iloc[j,0]])
		
#print(max_len_seq_table[17][41])
#print(similarity)
df = pd.DataFrame(similarity)
df.to_csv(output_file_name,header=None, index=None, sep='\t')




#similarity001 = [[] for i in range(max_id_of_people_001+1)]
#for i in range(len(input_seq[0])):
#	similarity001[input_seq.iloc[i,0]] = SimilarityOfLayer(input_seq,max_len_seq_table_001,idx_001,i)













# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 08:20:24 2018

@author: wth
"""

import pandas as pd
import numpy as np
import sys
import ast

input_file_name = sys.argv[1]
person_id = sys.argv[2]
mlength = sys.argv[3]

input_seq = pd.read_csv(input_file_name,header=None, encoding = "UTF-8", sep='\t')

max_id_of_people = max(input_seq[0])

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

	
def ExtendSequence(seq,i,j,tth):
#	print("seq:",seq)
#	print("enter ExtendSequence")
	res = []
	seqi = ast.literal_eval(input_seq.iloc[i,1])
	seqj = ast.literal_eval(input_seq.iloc[j,1])
	timei = ast.literal_eval(input_seq.iloc[i,3])
	timej = ast.literal_eval(input_seq.iloc[j,3])
	lastNode = seq[-1]
	for ii in range(lastNode.ipos+1,len(seqi)):
		deltai = sum(timei[lastNode.ipos:ii])
		for jj in range(lastNode.jpos+1,len(seqj)):
			deltaj = sum(timej[lastNode.jpos:jj])
#			print(seqi[ii],seqj[jj])
			if (seqi[ii]==seqj[jj]) & (abs(deltai-deltaj)<=tth):
				if(len(res)==0):
					res = [seq.copy()]
				res[-1].append(Node(seqi[ii],ii,jj))
#				print("s:",seq)
#				print("seq_extended:",res[-1])
#				print("s:",seq)
#				print()
				res.append(seq.copy())
					
	return res

def SequenceMatching(i,j,maxLength,tth):
	res = []
	step = 0
	seqi = ast.literal_eval(input_seq.iloc[i,1])
	seqj = ast.literal_eval(input_seq.iloc[j,1])
#	maxLength = max(len(seqi),len(seqj))
	
	#1_len_seq = []
	for x in range(len(seqi)):
		for y in range(len(seqj)):
			if (seqi[x]==seqj[y]):
				node = Node(seqi[x],x,y)
				res.append([node])
	
	step=1
#	print(res)
	while step<maxLength:
		extended_seqs = []
#		print("res:",res)
		for s in res:
#			print("s:",s,"len:",len(s))
			if len(s)==step:
				extended_seqs.append(ExtendSequence(s,i,j,tth))	

		
		for es in extended_seqs:
#			print("s:",s)
#			print("es:",es)
			res = res + es
#			print("res:",res)
		
		if len(res)>0:
			m = max(len(l) for l in res)
		else:
			break
#		print(m)
		step += 1
		res = [e for e in res if len(e)==m]
		
#		print("res:",res)
#	print("step:",step)
	return res

seqi = ast.literal_eval(input_seq.iloc[1,1])
seqj = ast.literal_eval(input_seq.iloc[1,1])

res = SequenceMatching(1,1,5,10800)
print("seqi:",seqi)
print("seqj:",seqj)
print("res:")
for r in res:
	print(r)

tth = 10800
max_length = int(mlength)
#max_len_seq_table = [[[] for i in range(max_id_of_people+1)] for i in range(max_id_of_people+1)]
seqi = ast.literal_eval(input_seq.iloc[int(person_id),1])
print("seq_"+str(input_seq.iloc[int(person_id),0])+":",seqi)
max_len_seq_table = [[] for i in range(max_id_of_people+1)]
for i in range(len(input_seq[0])):
	max_len_seq_table[input_seq.iloc[i,0]] = SequenceMatching(int(person_id),i,max_length,tth)
#	print(input_seq.iloc[int(person_id),0],"-",input_seq.iloc[i,0],":",max_len_seq_table[input_seq.iloc[i,0]])
#for i in range(len(input_seq[0])):
#	print(i)
#	for j in range(i,len(input_seq[0])):
#		max_len_seq_table[input_seq.iloc[i,0]][input_seq.iloc[j,0]] = SequenceMatching(i,j,max_length,tth)
#		print(input_seq.iloc[i,0],"-",input_seq.iloc[j,0],":",max_len_seq_table[input_seq.iloc[i,0]][input_seq.iloc[j,0]])
		
#print(max_len_seq_table[17][41])
df = pd.DataFrame(max_len_seq_table)
df.to_csv("person_"+str(person_id)+"_seq_table_"+str(max_length),header=None, index=None, sep='\t')

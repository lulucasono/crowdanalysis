# -*- coding: utf-8 -*-
"""
Created on Sun Nov 18 19:02:35 2018

@author: wth
"""

import pandas as pd
import numpy as np
import sys
import ast

input_file_name = sys.argv[1]
output_file_name = sys.argv[2]
person_id = sys.argv[3]
person_id2 = sys.argv[4]
max_length = sys.argv[5]
tth = sys.argv[6]
max_length = int(max_length)


input_seq = pd.read_csv(input_file_name,header=None, encoding = "UTF-8", sep='\t')
#input_seq = pd.read_csv("users_Seq_001",header=None, encoding = "UTF-8", sep='\t')
#input_seq2 = pd.read_csv("users_Seq_01",header=None, encoding = "UTF-8", sep='\t')
#input_seq3 = pd.read_csv("users_Seq_1",header=None, encoding = "UTF-8", sep='\t')

person_idx = input_seq.index[input_seq[0]==int(person_id)][0]
person_idx2 = input_seq.index[input_seq[0]==int(person_id2)][0]
#print(idx_001)
#idx_01 = input_seq2.index[input_seq2[0]==int(person_id)][0]
#print(idx_01)
#idx_1 = input_seq3.index[input_seq3[0]==int(person_id)][0]
#print(idx_1)

max_id_of_people = max(input_seq[0])
#max_id_of_people_01 = max(input_seq2[0])
#max_id_of_people_1 = max(input_seq3[0])

tth = int(tth)


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
#    print("Enter ExtendSequence")
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


#    print("res =", res)
#    input()
    return res

class LCS:
    def __init__(self, timei, timej, tth):
        self.timei = timei
        self.timej = timej
        self.tth = tth
        self.lcss = []
        self.state = [[0 for x in range(len(timej)+2)] for x in range(len(timei)+2)]
        
    def initState(self):
        self.state = [[0 for x in range(len(self.timej)+2)] for x in range(len(self.timei)+2)]
    
    def isLegal(self, seq):
        for x in range(1,len(seq)):
            deltai = sum(self.timei[seq[x].ipos:seq[x-1].ipos])
            deltaj = sum(self.timej[seq[x].jpos:seq[x-1].jpos])
            if abs(deltai-deltaj)>self.tth:
#                print(abs(deltai-deltaj))
                return False
        return True
    
    def getAllLcs(self, flag, seqi, lcs, maxSublen, i, j):
        while (i>0)&(j>0)&(len(lcs)<maxSublen):
#            print(i,j)
            if self.state[i][j]==0:
                self.state[i][j]=1
            else:
                break
            direction = flag[i][j]
            if direction == 2:
                if len(lcs)==0:
                    lcs.append(Node(seqi[i-1],i-1,j-1))
                    i -= 1
                    j -= 1  
                else:
                    lastNode = lcs[-1]
                    currentNode = Node(seqi[i-1],i-1,j-1)
                    deltai = sum(self.timei[currentNode.ipos:lastNode.ipos])
                    deltaj = sum(self.timej[currentNode.jpos:lastNode.jpos])
                    if abs(deltai-deltaj)<=self.tth:
                        lcs.append(Node(seqi[i-1],i-1,j-1))
                        i -= 1
                        j -= 1
                    else:
                        break
            else:
                if direction == 1:
                    i -= 1
                elif direction == 3:
                    j -= 1
                else:
#                    tmp = lcs.copy()
                    self.getAllLcs(flag, seqi, lcs.copy(), maxSublen, i-1, j)
#                    lcs = tmp
                    self.getAllLcs(flag, seqi, lcs.copy(), maxSublen, i, j-1)
#        if (len(lcs)==maxSublen) & self.isLegal(lcs):
#            self.lcss.append(lcs[::-1])
        if len(lcs)==maxSublen:
            self.lcss.append(lcs[::-1])
        del lcs




def SequenceMatching(users_seq, i, j ,tth):
    seqi = ast.literal_eval(users_seq.iloc[i,1])
    print(seqi)
    seqj = ast.literal_eval(users_seq.iloc[j,1])
    print(seqj)
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

    maxStep = min(max_length,max(max(c)))
    
#    print(pd.DataFrame(c))
#    print(maxStep)
#    print(pd.DataFrame(flag))

    ob_LCS = LCS(timei, timej, tth)
#    print(maxStep)
    a = np.asarray(c)
    sorted_values = (-a).argsort(axis=None, kind='mergesort')
    idx = np.unravel_index(sorted_values, a.shape)
    l = np.vstack(idx).T
    while (len(ob_LCS.lcss)==0) & (maxStep>0):
        ob_LCS.initState()
        print("======================")
        print(maxStep)
        n = 0
        for i in l:
            n += 1
            if c[i[0]][i[1]] >= maxStep:
                x = i[0]
                y = i[1]
                ob_LCS.getAllLcs(flag, seqi, [], maxStep, x, y)
            else:
                break
#        print(ob_LCS.lcss)
        print(n)
        print()
        maxStep -= 1
    
    print('max_length:',maxStep+1)            
#    print(ob_LCS.lcss)
    lcs_set = set(tuple(x) for x in ob_LCS.lcss)
    lcss = [ list(x) for x in lcs_set ]
#    print(lcss)
    return lcss

def SimilarityMeasureOfSequence(seq,i,j):
    seqNi = ast.literal_eval(input_seq.iloc[i,2])
    seqNj = ast.literal_eval(input_seq.iloc[j,2])
    if len(seq)==0:
        return 0

    sm = 0
    res = 0

    for node in seq:
        sm+=min(seqNi[node.ipos],seqNj[node.jpos])

    res = 2**(len(seq)-1) * sm

    return res

def SimilarityOfLayer(seq_list, i, j):
#    print("Enter SimilarityOfLayer")
    ni = input_seq.iloc[i,4]
    nj = input_seq.iloc[j,4]
    s = 0;
#    print(seq_list[input_seq.iloc[i,0]][input_seq.iloc[j,0]])
#    print(seq_list)
    for seq in seq_list[input_seq.iloc[i,0]][input_seq.iloc[j,0]]:
        s += SimilarityMeasureOfSequence(seq,i,j)
    s = s/(ni*nj)
    return s


max_len_seq_table = [[-1 for i in range(max_id_of_people+1)] for i in range(max_id_of_people+1)]
similarity = 0
print("=============================================")
#print("perseon id:",    pid)
max_len_seq_table[int(person_id)][int(person_id2)] = SequenceMatching(input_seq,person_idx,person_idx2,tth)
#print(max_len_seq_table[person_idx][person_idx2])
similarity = SimilarityOfLayer(max_len_seq_table, person_idx, person_idx2)
#print(person_id,"-",person_id2,":",pd.DataFrame(max_len_seq_table[int(person_id)][int(person_id2)]))
print(person_id,"-",person_id2,":",max_len_seq_table[int(person_id)][int(person_id2)])
print("Similarity("+str(person_id)+","+str(person_id2)+"):",similarity)

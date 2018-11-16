# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 18:37:42 2018

@author: wth
"""

import pandas as pd
import numpy as np
import sys
#import timeit

#start = timeit.default_timer()

input_file_name = sys.argv[1]
output_file_name = sys.argv[2]

input_table = pd.read_csv(input_file_name,header=None, encoding = "UTF-8", sep='\t')

# Remove noise
input_table = input_table[input_table[5] != -1]

number_of_sp = len(input_table[0])
print("number_of_sp:",number_of_sp)
number_of_people = len(input_table[0].drop_duplicates())
print("number_of_people:",number_of_people)

df = pd.DataFrame(index=range(number_of_people),columns=range(4))
df[1] = np.empty((len(df), 0)).tolist()
df[2] = np.empty((len(df), 0)).tolist()
df[3] = np.empty((len(df), 0)).tolist()

last_user_id = input_table.iloc[0,0]
last_cluster = input_table.iloc[0,5]
df.at[0,0] = input_table.iloc[0,0]
df.at[0,1].append(input_table.iloc[0,5])
df.at[0,2].append(1)
n=0

#print(df)

for i in range(number_of_sp):
	user_id = input_table.iloc[i,0]
	cluster = input_table.iloc[i,5]
	if(user_id==last_user_id):
		if(cluster==last_cluster):
			df.at[n,2][-1] += 1
		else:
			df.at[n,1].append(cluster)
			df.at[n,2].append(1)
			df.at[n,3].append(input_table.iloc[i,3]-input_table.iloc[i-1,4])
	else:
		n+=1
		df.at[n,0] = user_id
		df.at[n,1].append(input_table.iloc[i,5])
		df.at[n,2].append(1)
	last_user_id = user_id
	last_cluster = cluster	

	
df.to_csv(output_file_name,header=None, index=None, sep='\t')

#end = timeit.default_timer()
#print("Execution time: ",str(end-start),"s",sep='')
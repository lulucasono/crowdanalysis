# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 15:46:00 2018

@author: wth
"""

from __future__ import division
from math import sin, cos, sqrt, atan2, radians
#import numpy as np
import pandas as pd
import sys
from datetime import datetime

input_file_name = sys.argv[1]
output_file_name = sys.argv[2]

#Use utf-16 for file treated on Windows(usually encoded in UTF-16) 
input_table = pd.read_csv(input_file_name,header=None, encoding = "UTF-8", sep='\t')
#input_table = pd.read_csv(input_file_name,header=None, encoding = "UTF-16", sep='\t')

# Copy input to output
# output_table_1 = input_table.copy()
# output_table = pd.DataFrame(columns=['id','longitude', 'lattitude','enterTime','quitTime'])
dictionary = {}
rows_list = []

# Size of dataset
sLength = len(input_table[0])

#id user
id_user = input_table.iloc[0,0]

# Add a new column to mark the stay point cluster
#SP = pd.Series(-1 for i in range(sLength))
#output_table_1 = output_table_1.assign(SP=SP.values)

# Define the radius(kilometer) of the stay point
distThreh = 0.2;

# Define the time(second) of the stay point
timeThreh = 1800;

# Function to calcul distance between two points
def Distance(i,j):
	R = 6373.0
	lat1 = radians(input_table.iloc[i,3])
	lon1 = radians(input_table.iloc[i,2])
	lat2 = radians(input_table.iloc[j,3])
	lon2 = radians(input_table.iloc[j,2])
	dlon = lon2 - lon1
	dlat = lat2 - lat1	
	a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
	c = 2 * atan2(sqrt(a), sqrt(1 - a))	
	distance = R * c
	return distance

# Function to calcul average ixation
def MeanCood(i,j):
	sumLat = 0
	sumLong = 0
	for x in range(i,j+1):
		sumLat += input_table.iloc[x,3]
		sumLong += input_table.iloc[x,2]
	return sumLong/(j+1-i),sumLat/(j+1-i)

# Main loop
i=0
n=0
while i<sLength:
	#print(i)
	j=i+1
	Token=0
	while j<sLength:
		dist = Distance(i,j)
		if dist>distThreh:
			timej = input_table.iloc[j,1]
			timei = input_table.iloc[i,1]
			timestampj = 0
			timestampi = 0
			if len(timej)>19:
				timestampj = datetime.strptime(timej, "%Y-%m-%d %H:%M:%S.%f").timestamp()
			else:
				timestampj = datetime.strptime(timej, "%Y-%m-%d %H:%M:%S").timestamp()
			if len(timei)>19:
				timestampi = datetime.strptime(timei, "%Y-%m-%d %H:%M:%S.%f").timestamp()
			else:
				timestampi = datetime.strptime(timei, "%Y-%m-%d %H:%M:%S").timestamp()
			deltaT = timestampj - timestampi
			if deltaT>timeThreh:
				#print(n)
				Token2 = 0
				a,b = MeanCood(i,j)
				#print(a,b,output_table_1.ix[i,'timestamp'],output_table_1.ix[j,'timestamp'])                
				#output_table_2 = output_table_2.append([{'index':n,'longitude':a,'lattitude':b,'enterTime':output_table_1.ix[i,'timestamp']-output_table_1.ix[i,'timestamp'],'count':1}], ignore_index=True)
				dictionary[n] = [id_user,a,b,timestampi,timestampj]
				n+=1
				#print(dictionary)
				i = j
				Token = 1
			break
		j+=1
	if Token!=1:
		i+=1

#print(output_table_1.to_string())
#print(output_table_2.to_string())
#output_table_1.to_csv('./out1.csv')
output_table = pd.DataFrame.from_dict(dictionary,  orient='index', columns=['id','longitude', 'lattitude','enterTime','quitTime']);
output_table.to_csv(output_file_name,header=None, index=None, sep='\t')

	

# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 15:06:39 2018

@author: wth
"""

import pandas as pd
import numpy as np
import sys
from sklearn.cluster import DBSCAN

input_file_name = sys.argv[1]
output_file_name = sys.argv[2]
min_distance = float(sys.argv[3])

input_table = pd.read_csv(input_file_name,header=None, encoding = "UTF-8", sep='\t')

data = input_table.iloc[:, 1:3]

kms_per_radian = 6371.0088
epsilon = min_distance / kms_per_radian
db = DBSCAN(eps=epsilon, min_samples=2, algorithm='ball_tree', metric='haversine').fit(np.radians(data))
input_table[len(input_table.columns)] = db.labels_
print (max(db.labels_))

input_table.to_csv(output_file_name,header=None, index=None, sep='\t')
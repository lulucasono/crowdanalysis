# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 23:32:18 2018

@author: wth
"""

from __future__ import division
import pandas as pd
import sys

input_file_name = sys.argv[1]
output_file_name = sys.argv[2]

input_table = pd.read_csv(input_file_name,header=None, encoding = "UTF-8", sep='\t')

output = input_table.iloc[::20, :]
#print(output)
output.to_csv(output_file_name,header=None, index=None, sep='\t')

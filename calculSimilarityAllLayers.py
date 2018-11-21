import pandas as pd
#import numpy as np
import sys
import ast

output_file_name='similarityAllLayers'

input_file_name='similarity_1'
input_seq_1 = pd.read_csv(input_file_name,header=None, encoding = "UTF-8", sep='\t')
input_file_name='similarity_01'
input_seq_01 = pd.read_csv(input_file_name,header=None, encoding = "UTF-8", sep='\t')
input_file_name='similarity_001'
input_seq_001 = pd.read_csv(input_file_name,header=None, encoding = "UTF-8", sep='\t')


dfsum=input_seq_1+2*input_seq_01+4*input_seq_001
dfsum.to_csv(output_file_name,header=None, index=None, sep='\t')

import pandas as pd
import sys

n = int(sys.argv[1])
input_table = ['' for i in range(n)]
for i in range(n):
    input_table[i] = sys.argv[2+i]
output_file_name=sys.argv[2+n]

dfs = []

ssum = 0
for i in range(n):
    dfs.append(pd.read_csv(input_table[i],header=None, encoding = "UTF-8", sep='\t'))

dfsum=dfs[i]
for i in range(1,n):
    dfsum += 2**i*dfs[i]
    
dfsum.to_csv(output_file_name,header=None, index=None, sep='\t')

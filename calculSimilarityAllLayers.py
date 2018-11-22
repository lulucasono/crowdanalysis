import pandas as pd
import sys

input1 = sys.argv[1]
input2 = sys.argv[2]
input3 = sys.argv[3]
output_file_name=sys.argv[4]

df1 = pd.read_csv(input1,header=None, encoding = "UTF-8", sep='\t')
df2 = pd.read_csv(input2,header=None, encoding = "UTF-8", sep='\t')
df3 = pd.read_csv(input3,header=None, encoding = "UTF-8", sep='\t')

dfsum=df1+2*df2+4*df3
dfsum.to_csv(output_file_name,header=None, index=None, sep='\t')

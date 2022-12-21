import pandas as pd

df = pd.read_csv("bold_Netherlands.fasta", sep=',', header=0)
print(df)
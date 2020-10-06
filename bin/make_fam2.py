#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys

csv_file = sys.argv[1]
df = pd.read_csv(csv_file, sep=',', header=0)

# remove vcf column
df = df.drop(['vcf'], axis=1)
# make Family ID column
df['FID'] = np.arange(len(df)) + 1
# rename vcf with IID
df = df.rename({'sampleid': 'IID'}, axis=1)
# make paternal & maternal ID columns
df['PAT'] = 0
df['MAT'] = 0

# reorder columns
cols = list(df.columns.values)
new_cols = []
main_cols = ['FID','IID','PAT','MAT','SEX']
for col in main_cols:
    cols.remove(col)
    new_cols.append(col)

for col in cols:
    new_cols.append(col)

df = df[new_cols]

df.to_csv('sample.phe', sep=' ', index=False)

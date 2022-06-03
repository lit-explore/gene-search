#!/bin/env python
"""
Use Biothings API to construct a gene symbol -> entrez gene id mapping
kh (2022/03/23)
"""
import time
import pandas as pd
from biothings_client import get_client

mg = get_client("gene")

# load list of hugo gene symbols
# https://www.genenames.org/download/custom/
hugo = pd.read_csv("/data/ref/hugo/genenames_2020-08-08.tsv", sep='\t')

# split gene symbols into batches
# https://stackoverflow.com/a/312464/554531
def chunks(lst, n):
  for i in range(0, len(lst), n):
      yield lst[i:i + n]

batches = list(chunks(hugo['Approved symbol'].values, 100))

res = []

for i, batch in enumerate(batches):
    print(f"Querying Biothings ({i + 1}/{len(batches)})")
    mapping = mg.querymany(batch, scopes="symbol", 
                           species=9606, fields="entrezgene", as_dataframe=True)

    mapping = mapping[~mapping.entrezgene.isna()]

    res.append(mapping)

# drop all columns except the index (symbol) and entrez ids
df = pd.concat([x['entrezgene'] for x in res]).to_frame()

df = df.reset_index().rename(columns={"query": "symbol"})

df.to_csv("../api/app/data/symbol_entrez.tsv", sep='\t', index=False)

#!/bin/env python
"""
Use Biothings API to construct a ensgene -> entrez gene id mapping
kh (2022/06/02)
"""
import time
import pandas as pd
from biothings_client import get_client

mg = get_client("gene")

# load list of ensembl gene identifiers
# https://www.genenames.org/download/custom/
with open("/data/ref/hugo/ensgenes_2022-06-02.txt") as fp:
    ensgenes = [x.strip() for x in fp.readlines()]

# split query into batches
# https://stackoverflow.com/a/312464/554531
def chunks(lst, n):
  for i in range(0, len(lst), n):
      yield lst[i:i + n]

batches = list(chunks(ensgenes, 100))

res = []

for i, batch in enumerate(batches):
    print(f"Querying Biothings ({i + 1}/{len(batches)})")
    mapping = mg.querymany(batch, 
                           scopes="ensembl.gene", 
                           species=9606, fields="entrezgene", 
                           as_dataframe=True)

    mapping = mapping[~mapping.entrezgene.isna()]

    res.append(mapping)

# drop all columns except the index (ensgene) and entrez ids
df = pd.concat([x['entrezgene'] for x in res]).to_frame()

df = df.reset_index().rename(columns={"query": "ensgene"})

# drop duplicated ensgene entries;
# only a single one found in testing: ENSG00000286105
df = df.drop_duplicates('ensgene')

df.to_csv("../api/app/data/ensgene_entrez.tsv", sep='\t', index=False)

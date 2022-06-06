#!/bin/env python
"""
Pubtator x Gene Search
KH (Oct2021)

- [ ] instead of evaluating all pubmed articles with >= 1 gene match,
      choose a cutoff that scales with the # input genes; this should help speed things
      up significantly for larger input queries
- [ ] cluster summary (genes assoc. with each cluster?)
- [ ] cache pmid citation lookups to avoid re-querying server?..
    - [ ] store mapping from gene queries to files
- [ ] add k-means "k" param to api?

"""
import datetime
import json
import os
import sys
import logging
import numpy as np
import pandas as pd
from collections import OrderedDict
from biothings_client import get_client
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from typing import Optional, List

from genesearch.pubmed import query_article_info, build_article_gene_comat 
from genesearch.clustering import cluster_articles
from genesearch.network import build_network

app = FastAPI()

# setup logger
logging.basicConfig(stream=sys.stdout, format='%(asctime)s [%(levelname)s] %(message)s', 
                    level=logging.INFO)
logger = logging.getLogger('gene-search')

# CORS
origins = [
    "http://localhost:81",
    "http://lit.biodat.io"
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# minimum p-value cutoff to use when compute weighted article scores
# TODO: convert to api parameter..
MIN_PVAL = 1e-100

# load gene symbol -> entrez id mapping (created: march 23, 2022)
# n = 43,491 gene symbols
#            entrezgene
#  symbol
#  A1BG               1
#  A1BG-AS1      503538
#  A1CF           29974
symbol_entrez_mapping = pd.read_csv("data/symbol_entrez.tsv", 
                                    sep='\t').set_index('symbol')

# load ensembl gene id -> entrez id mapping (created: june 02, 2022)
ensgene_entrez_mapping = pd.read_csv("data/ensgene_entrez.tsv", 
                                     sep='\t').set_index('ensgene')

# load disease mesh term -> pmid mapping
# n ~ 11,000 MeSH terms
disease_pmid_infile = "/data/pmids/disease-pmids.json"

with open(disease_pmid_infile) as fp:
    mesh_pmids = json.loads(fp.read())

# load entrez id -> pmid mapping
# generated from pubtator data
# n ~ 26,300 entrez ids
entrez_pmid_infile = "/data/pmids/gene-pmids.json"

with open(entrez_pmid_infile) as fp:
    entrez_article_mapping = json.loads(fp.read())

# load pmid -> unique gene counts
# number of unique gene mentions for each pubmed article
# n ~ 4.35m pubmed ids
gene_counts_infile = "/data/gene-counts/pmid_gene_counts.feather"
gene_counts = pd.read_feather(gene_counts_infile).set_index("pmid")

# create output dir, if it doesn't already exist
if not os.path.exists("/data/gene-search/"):
    os.makedirs("/data/gene-search/", mode=0o755)

@app.post("/api/")
async def query(genes: str, keyType: str, pvalues: Optional[str] = None, 
                disease: Optional[str] = None, 
                max_articles: Optional[int] = 100):
    # split list of genes
    input_genes = genes.split(",")

    # validate key type
    if keyType not in ['symbol', 'ensgene']:
        return ({"error": "Invalid keyType specified!"})

    # validate disease MeSH term, if specified
    if disease is not None:
        if disease not in mesh_pmids.keys():
            return({"error": "Unable to find specified disease MeSH ID"})

    # validate p-value input length, and create series for easier masking later on
    if pvalues is None:
        input_pvalues = [None for i in range(len(input_genes))]
    else:
        input_pvalues = pvalues.split(",")

        if len(input_pvalues) != len(input_genes):
            return({"error": "Gene symbol and P-value inputs have different lengths!"})

        # convert pvalues to numeric
        input_pvalues = [np.float64(pval) for pval in input_pvalues]

    # combine input genes and p-values into a single dataframe, in order to make them
    # easier to transform together;
    # if no p-values specified, a column of all missing values will be used as a
    # placeholder
    query_df = pd.DataFrame({
        "gene": input_genes,
        "pvalue": input_pvalues
    })

    # create a dict mapping from gene symbol -> p-value
    #pvalues_series = pd.Series([np.float64(pval) for pval in pvalues.split(",")])

    # if disease specified, limit to articles pertaining to that disease
    pmids_allowed = None

    if disease is not None:
        pmids_allowed = mesh_pmids[disease]

    # if ensembl gene ids or symbols provided, convert to entrez ids
    if keyType in ['ensgene', 'symbol']:
        # determine gene mapping to use
        if keyType == 'ensgene':
            gene_mapping = ensgene_entrez_mapping
        else:
            gene_mapping = symbol_entrez_mapping

        # exclude any genes that can't be mapped to entrez ids
        mask = query_df.gene.isin(gene_mapping.index)
        query_df = query_df[mask]

        # add entrez ids column
        entrez_ids = [str(x) for x in gene_mapping.loc[query_df.gene].entrezgene.values]
        query_df.insert(0, "entrezgene", entrez_ids)
    else:
        query_df.insert(0, "entrezgene", input_genes)
    
    # exclude any genes that aren't present in the pubtator data
    valid_entrez_ids = entrez_article_mapping.keys()

    mask = query_df.entrezgene.isin(valid_entrez_ids)
    query_df = query_df[mask]

    # to speed things up, get subset of gene/article mapping dict for genes of interest
    entrez_article_mapping_subset = {}

    entrez_ids = query_df.entrezgene.values

    for entrez_id in entrez_ids:
        entrez_article_mapping_subset[entrez_id] = entrez_article_mapping[entrez_id]
        entrez_article_mapping_subset[entrez_id]

    # get binary <article x gene> matrix
    comat = build_article_gene_comat(entrez_ids, entrez_article_mapping_subset, 
                                     pmids_allowed)

    # devel: store pmid x gene mat..
    #  now = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")
    #  comat.reset_index().rename(columns={"index": "pmid"}).to_feather(f"/data/gene-search/{now}.feather")

    # get number of _target_ genes matched in each article
    num_matched = comat.sum(axis=1)

    # get the _total_ number of unique genes associated with each article
    article_total_genes = gene_counts.loc[num_matched.index].num

    # ratio of matched genes to total genes
    ratio_matched = num_matched / article_total_genes

    # if p-values were provided, weight gene contributions by -log10 (p-values)
    if pvalues is not None:
        dat_weighted = comat.copy()

        # create a dict mapping from gene symbol -> p-value
        #  target_pvals = {}
        #
        #  for i, pval in enumerate(query_df.pvalue.values):
        #      target_pvals[target_symbols[i]] = np.float64(pval)

        # sanity check..
        if not np.all(query_df.entrezgene == comat.columns):
            raise Exception("ID mismatch during weighted co-mat generation!")

        # multiply columns of co-ocurrence matrix by p-value based weights
        for i, gene in query_df.iterrows():
            dat_weighted.iloc[:, i] *= -np.log10(np.max([MIN_PVAL, gene.pvalue]))

        # rescale to [0, max(num_matched)]
        max_matches = num_matched.max().max()
        max_weight = dat_weighted.max().max()

        dat_weighted = (dat_weighted / max_weight) * max_matches

        # get total scores for each article
        article_scores = dat_weighted.sum(axis=1)
    else:
        article_scores = num_matched

    #
    # score = <matched>^2 * <matched/total>
    #
    # this way, articles mentioning a large number of target genes are prioritized,
    # while penalizing articles which simply include a large number of genes overall
    #
    article_scores = article_scores** 2 * ratio_matched

    article_scores = article_scores.sort_values(ascending=False).head(max_articles)

    # get subset of article x gene matrix corresponding to the top N hits and cluster
    # articles based on similarity in the genes mentioned..
    comat = comat[comat.index.isin(article_scores.index)]

    # drop genes which are not present in any articles, after filtering..
    mask = comat.sum() > 0

    if (~mask).sum() > 0:
        dropped_genes = ", ".join(sorted(comat.columns[~mask]))

        logger.info(
            f"Dropping genes that don't appear in any of the top-scoring articles: {dropped_genes}"
        )
        comat = comat.loc[:, mask]
    else:
        dropped_genes = ""

    logger.info(f"Clustering {comat.shape[0]} articles...")
    clusters = cluster_articles(comat)

    # generate network
    logger.info("Building article network..")
    net = build_network(comat)

    # query pubmed api for article info
    citations = query_article_info(list(comat.index))

    # cluster color map
    cmap = [
        "#aef980",
        "#4ff0b9",
        "#27dde7",
        "#83c0f3",
        "#c79cd8",
        "#cb7ffa",
        "#ea3e8b",
        "#e67da1",
        "#fa8e7f",
        "#f2e283",
    ]

    # generate section with cluster information
    clust_info = {}

    logger.info("Finalizing results...")

    # swap indices in gene id mapping to make it easier to convert back to the original
    # key type
    gene_mapping = gene_mapping.reset_index().set_index('entrezgene')

    for clust in set(clusters.values()):
        # get pmids associated with cluster
        clust_dict = dict(filter(lambda x: x[1] == clust, clusters.items()))
        clust_pmids = clust_dict.keys()

        # get submat associated with cluster, and count genes
        clust_comat = comat[comat.index.isin(clust_pmids)]

        # get within-cluster target gene counts
        clust_counts = clust_comat.sum()

        # if original query used gene symbols / ensgenes, convert back
        if keyType in ['ensgene', 'symbol']:
            orig_ids = gene_mapping.loc[clust_counts.index.astype('int')][keyType].values
            clust_counts.index = orig_ids

        count_dict = clust_counts.sort_values(ascending=False).to_dict(into=OrderedDict)

        clust_info[clust] = {
            "color": cmap[clust],
            "genes": count_dict,
            "num_articles": len(clust_pmids),
        }

    scores_dict = article_scores.to_dict()

    dat_subset_sums = num_matched.loc[comat.index]

    # generate result json, including gene ids, etc. for each article
    articles = []

    for pmid, num_genes_matched in dat_subset_sums.items():

        # retrieve the associated gene names
        pmid_genes = list(comat.loc[pmid][comat.loc[pmid] == 1].index)

        # convert back to original key type, if needed
        if keyType in ['ensgene', 'symbol']:
            pmid_genes = gene_mapping.loc[[int(x) for x in pmid_genes]][keyType].values

        articles.append(
            {
                "id": pmid,
                "genes": ", ".join(pmid_genes),
                "num_matched": num_genes_matched,
                "num_total": gene_counts.loc[pmid].num.item(),
                "score": scores_dict[pmid],
                "citation": citations[pmid],
                "cluster": clusters[pmid],
                "color": cmap[clusters[pmid]],
            }
        )

    # sort articles by score and add to "network" object
    articles.sort(key=lambda x: x["score"], reverse=True)
    net["nodes"] = articles

    info = {"not_present_in_top_articles": dropped_genes}

    query_dict = {"genes": input_genes}

    if pvalues is not None:
        query_dict["pvalues"] = input_pvalues 

    res = {"query": query_dict, "clusters": clust_info, "network": net, "info": info}

    # testing..
    now = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")

    with open(f"/data/gene-search/{now}.json", "w") as fp:
        json.dump(res, fp)

    print("Done!")

    return res


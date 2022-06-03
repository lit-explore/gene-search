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

from genesearch.pubmed import query_article_info, get_pmid_symbol_comat
from genesearch.clustering import cluster_articles
from genesearch.network import build_network

app = FastAPI()
#app = FastAPI(root_path="/api")

# setup logger
logging.basicConfig(stream=sys.stdout, format='%(asctime) [%(levelname)s] %(message)s', 
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
gene_mapping_infile = "data/symbol_entrez.tsv"

with open(gene_mapping_infile) as fp:
    gene_mapping = pd.read_csv(gene_mapping_infile, sep='\t').set_index('symbol')


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
    entrez_pmids = json.loads(fp.read())

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

    # validate p-value input length, if specified
    if pvalues is not None:
        if len(pvalues.split(",")) != len(input_genes):
            return({"error": "Gene symbol and P-value inputs have different lengths!"})

    # if disease specified, limit to articles pertaining to that disease
    pmids_allowed = None

    if disease is not None:
        pmids_allowed = mesh_pmids[disease]

    # get binary <pmid x gene> matrix
    pmid_symbol_comat = get_pmid_symbol_comat(input_genes, gene_mapping,
                                              entrez_pmids, pmids_allowed)

    # devel: store pmid x gene mat..
    now = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")
    pmid_symbol_comat.reset_index().rename(columns={"index": "pmid"}).to_feather(f"/data/gene-search/{now}.feather")

    # get number of genes matched in each article
    num_matched = pmid_symbol_comat.sum(axis=1)

    # get the _total_ number of genes associated with each article
    pmid_gene_counts = gene_counts.loc[num_matched.index].num

    # ratio of matched genes to total genes
    ratio_matched = num_matched / pmid_gene_counts

    # if p-values were provided, weight gene contributions by -log10 (p-values)
    if pvalues is not None:
        dat_weighted = pmid_symbol_comat.copy()

        # create a dict mapping from gene symbol -> p-value
        target_pvals = {}

        for i, pval in enumerate(pvalues.split(",")):
            target_pvals[target_symbols[i]] = np.float64(pval)

        # multiply columns of co-ocurrence matrix by p-value based weights
        for i, symbol in enumerate(dat_weighted.columns):
            dat_weighted.iloc[:, i] *= -np.log10(np.max([MIN_PVAL, target_pvals[symbol]]))

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
    pmid_symbol_comat_top = pmid_symbol_comat[pmid_symbol_comat.index.isin(article_scores.index)]

    # drop genes which are not present in any articles, after filtering..
    mask = pmid_symbol_comat_top.sum() > 0

    dropped_genes = ", ".join(sorted(pmid_symbol_comat_top.columns[~mask]))

    logger.info(
        f"Dropping genes that don't appear in any of the top-scoring articles: {dropped_genes}"
    )
    pmid_symbol_comat_top = pmid_symbol_comat_top.loc[:, mask]

    logger.info(f"Clustering {pmid_symbol_comat_top.shape[0]} articles...")
    clusters = cluster_articles(pmid_symbol_comat_top)

    # generate network
    logger.info("Building article network..")
    net = build_network(pmid_symbol_comat_top)

    # query pubmed api for article info
    citations = query_article_info(list(pmid_symbol_comat_top.index))

    # cluster color map
    # todo: add check & interpolate colors if more are needed..
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

    for clust in set(clusters.values()):
        # get pmids associated with cluster
        clust_dict = dict(filter(lambda x: x[1] == clust, clusters.items()))
        clust_pmids = clust_dict.keys()

        # get submat associated with cluster, and count genes
        clust_submat = pmid_symbol_comat_top[pmid_symbol_comat_top.index.isin(clust_pmids)]

        clust_counts = clust_submat.sum()
        count_dict = clust_counts.sort_values(ascending=False).to_dict(into=OrderedDict)

        clust_info[clust] = {
            "color": cmap[clust],
            "genes": count_dict,
            "num_articles": len(clust_pmids),
        }

    scores_dict = article_scores.to_dict()

    dat_subset_sums = num_matched.loc[pmid_symbol_comat_top.index]

    # generate result json, including gene ids, etc. for each article
    articles = []

    for pmid, num_genes in dat_subset_sums.items():

        # retrieve the associated gene names
        pmid_genes = list(pmid_symbol_comat_top.loc[pmid][pmid_symbol_comat_top.loc[pmid] == 1].index)

        articles.append(
            {
                "id": pmid,
                "genes": ", ".join(pmid_genes),
                "num_matched": num_genes,
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

    query_dict = {"genes": target_symbols}

    if pvalues is not None:
        query_dict["pvalues"] = [np.float64(pval) for pval in pvalues.split(",")]

    res = {"query": query_dict, "clusters": clust_info, "network": net, "info": info}

    # testing..
    with open(f"/data/gene-search/{now}.json", "w") as fp:
        json.dump(res, fp)

    print("Done!")

    return res


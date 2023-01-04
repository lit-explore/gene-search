"""
gene-search
KH (Oct2021)
"""
import datetime
import json
import os
import sys
import logging
from collections import OrderedDict
from typing import Optional
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
import numpy as np
import pandas as pd

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
MIN_PVAL = 1e-100

# cluster color map
CMAP = [
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
with open("/data/pmids/disease-pmids.json", "rt", encoding="utf-8") as fp:
    mesh_pmids = json.loads(fp.read())

# load entrez id -> pmid mapping
# generated from pubtator data
# n ~ 26,300 entrez ids
with open("/data/pmids/gene-pmids.json", "rt", encoding="utf-8") as fp:
    entrez_article_mapping = json.loads(fp.read())

# load pmid -> unique gene counts
# number of unique gene mentions for each pubmed article
# n ~ 4.35m pubmed ids
gene_counts = pd.read_feather("/data/gene-counts/pmid_gene_counts.feather").set_index("pmid")

# create output dir, if it doesn't already exist
if not os.path.exists("/data/gene-search/"):
    os.makedirs("/data/gene-search/", mode=0o755)

def create_cluster_info(clusters:dict, comat:pd.DataFrame, key_type:str, 
                        gene_mapping:pd.DataFrame|None):
    """
    Creates cluster info section of result
    """
    # generate section with cluster information
    clust_info = {}

    for clust in set(clusters.values()):
        # get pmids associated with cluster
        clust_dict = dict(filter(lambda x: x[1] == clust, clusters.items()))
        clust_pmids = clust_dict.keys()

        # get submat associated with cluster, and count genes
        clust_comat = comat[comat.index.isin(clust_pmids)]

        # get within-cluster target gene counts
        clust_counts = clust_comat.sum()

        # if original query used gene symbols / ensgenes, convert back
        if gene_mapping is not None:
            ind = clust_counts.index.astype('int')
            clust_counts.index = gene_mapping.loc[ind][key_type].values

        count_dict = clust_counts.sort_values(ascending=False).to_dict(into=OrderedDict)

        clust_info[clust] = {
            "color": CMAP[clust],
            "genes": count_dict,
            "num_articles": len(clust_pmids),
        }

@app.post("/api/")
async def query(genes: str, key_type: str, pvalues: Optional[str] = None,
                disease: Optional[str] = None,
                max_articles: Optional[int] = 100):
    # validate key type
    if key_type not in ['symbol', 'ensgene']:
        return {"error": "Invalid key_type specified!"}

    # validate disease MeSH term, if specified
    if disease is not None:
        if disease not in mesh_pmids.keys():
            return {"error": "Unable to find specified disease MeSH ID"}

    # split list of genes
    input_genes = genes.split(",")

    # create dict with query information to be included in result
    query_dict = {"genes": input_genes}

    # validate p-value input length, and create series for easier masking later on
    if pvalues is None:
        input_pvalues = [None for i in range(len(input_genes))]
        query_dict["pvalues"] = input_pvalues
    else:
        input_pvalues = pvalues.split(",")

        if len(input_pvalues) != len(input_genes):
            return {"error": "Gene symbol and P-value inputs have different lengths!"}

        # convert pvalues to numeric
        input_pvalues = [np.float64(pval) for pval in input_pvalues]

    # combine input genes and p-values into a single dataframe, in order to make them
    # easier to transform together;
    # if no p-values specified, a column of all missing values will be used as a
    # placeholder
    input_df = pd.DataFrame({
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
    gene_mapping = None

    if key_type in ['ensgene', 'symbol']:
        # determine gene mapping to use
        if key_type == 'ensgene':
            gene_mapping = ensgene_entrez_mapping
        else:
            gene_mapping = symbol_entrez_mapping

        # exclude any genes that can't be mapped to entrez ids
        mask = input_df.gene.isin(gene_mapping.index)
        input_df = input_df[mask]

        # add entrez ids column
        entrez_ids = [str(x) for x in gene_mapping.loc[input_df.gene].entrezgene.values]
        input_df.insert(0, "entrezgene", entrez_ids)
    else:
        input_df.insert(0, "entrezgene", input_genes)

    # exclude any genes that aren't present in the pubtator data
    valid_entrez_ids = entrez_article_mapping.keys()

    mask = input_df.entrezgene.isin(valid_entrez_ids)
    input_df = input_df[mask]

    # to speed things up, get subset of gene/article mapping dict for genes of interest
    entrez_article_mapping_subset = {}

    entrez_ids = input_df.entrezgene.values

    for entrez_id in entrez_ids:
        entrez_article_mapping_subset[entrez_id] = entrez_article_mapping[entrez_id]
        entrez_article_mapping_subset[entrez_id]

    # get binary <article x gene> matrix
    comat = build_article_gene_comat(entrez_ids, entrez_article_mapping_subset,
                                     pmids_allowed)

    # devel: store pmid x gene mat..
    #  now = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")
    #  comat.reset_index().rename(columns={"index": "pmid"})
    #    .to_feather(f"/data/gene-search/{now}.feather")

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
        #  for i, pval in enumerate(input_df.pvalue.values):
        #      target_pvals[target_symbols[i]] = np.float64(pval)

        # sanity check..
        if not np.all(input_df.entrezgene == comat.columns):
            raise Exception("ID mismatch during weighted co-mat generation!")

        # multiply columns of co-ocurrence matrix by p-value based weights
        for i, gene in input_df.iterrows():
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

    logger.info("Finalizing results...")

    # swap indices in gene id mapping to make it easier to convert back to the original
    # key type
    if gene_mapping is not None:
        gene_mapping = gene_mapping.reset_index().set_index('entrezgene')

    # create cluster info section
    clust_info = create_cluster_info(clusters, comat, key_type, gene_mapping)

    scores_dict = article_scores.to_dict()

    dat_subset_sums = num_matched.loc[comat.index]

    # generate result json, including gene ids, etc. for each article
    articles = []

    for pmid, num_genes_matched in dat_subset_sums.items():

        # retrieve the associated gene names
        pmid_genes = list(comat.loc[pmid][comat.loc[pmid] == 1].index)

        # convert back to original key type, if needed
        if key_type in ['ensgene', 'symbol']:
            pmid_genes = gene_mapping.loc[[int(x) for x in pmid_genes]][key_type].values

        articles.append(
            {
                "id": pmid,
                "genes": ", ".join(pmid_genes),
                "num_matched": num_genes_matched,
                "num_total": gene_counts.loc[pmid].num.item(),
                "score": scores_dict[pmid],
                "citation": citations[pmid],
                "cluster": clusters[pmid],
                "color": CMAP[clusters[pmid]],
            }
        )

    # sort articles by score and add to "network" object
    articles.sort(key=lambda x: x["score"], reverse=True)
    net["nodes"] = articles

    info = {"not_present_in_top_articles": dropped_genes}

    res = {"query": query_dict, "clusters": clust_info, "network": net, "info": info}

    # store query results as json
    now = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")

    outfile = f"/data/gene-search/{now}.json"

    with open(outfile, "w", encoding="utf-8") as fp:
        json.dump(res, fp)

    # store query results as tsv
    pd.DataFrame(articles).to_csv(outfile.replace('.json', '.tsv'), sep='\t')

    return res

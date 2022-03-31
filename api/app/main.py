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
import numpy as np
import pandas as pd
import urllib.request
from collections import OrderedDict
from biothings_client import get_client
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from sklearn.cluster import KMeans
from typing import Optional

app = FastAPI()

# CORS
origins = [
    "http://localhost:81",
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

# load gene symbol -> entrez id mapping (march 23, 2022)
# n = 43,491 gene symbols
gene_mapping_infile = "data/symbol_entrez.tsv"

with open(gene_mapping_infile) as fp:
    gene_mapping = pd.read_csv(gene_mapping_infile, sep='\t').set_index('symbol')

# load disease mesh term -> pmid mapping
# n ~ 11,000 MeSH terms
disease_pmid_infile = "/data/pmids/disease-pmids.json"

with open(disease_pmid_infile) as fp:
    mesh_pmids = json.loads(fp.read())

# load entrez -> pmid mapping
# n ~ 26,300 entrez ids
entrez_pmid_infile = "/data/pmids/gene-pmids.json"

with open(entrez_pmid_infile) as fp:
    entrez_pmids = json.loads(fp.read())

# pmid -> total # unique gene mentions for each article
# n ~ 4.35m pubmed ids
gene_counts_infile = "/data/gene-counts/pmid_gene_counts.feather"

gene_counts = pd.read_feather(gene_counts_infile).set_index("pmid")

# create output dir, if it doesn't already exist
if not os.path.exists("/data/gene-search/"):
    os.makedirs("/data/gene-search/", mode=0o755)

def get_pmid_symbol_comat(target_symbols, disease=None):
    """Generates an <article x gene> matrix for genes of interest"""
    # check for any genes that could not be mapped and remove them..
    target_symbols = pd.Series(target_symbols)

    # exclude genes that could not be mapped, or whose mapped entrez ids are
    # not present in the pubtator dataset
    mask = target_symbols.isin(gene_mapping.index)
    target_symbols = target_symbols[mask]

    # convert gene target_symbols to entrez ids (entrez_pmids is indexed by string ids)
    target_entrez = gene_mapping.loc[target_symbols].entrezgene.values
    target_entrez = [str(x) for x in target_entrez]

    # exlude genes for which no entrez is present in the entrez -> pmid mapping
    mask = pd.Series([x in entrez_pmids.keys() for x in target_entrez])
    target_entrez = pd.Series(target_entrez)[mask].values

    target_symbols = target_symbols[mask.values]

    # subset pmids to get only genes of interest
    entrez_pmids_subset = {}

    for entrez_id in target_entrez:
        entrez_pmids_subset[entrez_id] = entrez_pmids[entrez_id]
        entrez_pmids_subset[entrez_id]

    # create a list of lists containing all pmids associated with >= target gene
    pmid_lists = [entrez_pmids_subset[entrez_id] for entrez_id in entrez_pmids_subset]

    # flatten, sort, remove duplicates, and convert to a series
    matching_pmids = pd.Series(sorted(set(sum(pmid_lists, []))))

    num_matches = len(matching_pmids)
    print(f"Found {num_matches} articles with one or more of the genes of interest..")

    # if disease specified, filter to include only articles relating to that disease
    if disease is not None:
        num_before = num_matches
        matching_pmids = matching_pmids[matching_pmids.isin(mesh_pmids[disease])]
        num_after = len(matching_pmids)

        print(f"Excluding {num_before - num_after}/{num_before} articles not matching specified disease.")

    # convert to a binary <pmid x gene> co-occurrence matrix
    pmid_symbol_rows = []

    for i, gene in enumerate(target_entrez):
        row = matching_pmids.isin(entrez_pmids_subset[gene]).astype(np.int64)
        pmid_symbol_rows.append(row)

    pmid_symbol_mat = pd.concat(pmid_symbol_rows, axis=1)

    # convert boolean to numeric (0|1)
    #  pmid_symbol_mat.replace({False: 0, True: 1}, inplace=True)

    pmid_symbol_mat.columns = target_symbols
    pmid_symbol_mat.index = matching_pmids

    return pmid_symbol_mat

# helper function chunk a list
# https://stackoverflow.com/a/312464/554531
def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i : i + n]

def query_article_info(pmids):
    """Uses PubMed API to retrieve additional information for each article hit"""
    base_url = "https://api.ncbi.nlm.nih.gov/lit/ctxp/v1/pubmed/?format=citation&id="

    citations = {}

    pmid_chunks = list(chunks(pmids, 50))

    for i, chunk in enumerate(pmid_chunks):
        print(f"Querying PubMed API for citation info.. ({i + 1}/{len(pmid_chunks)})")

        url = base_url + ",".join([str(x) for x in chunk])

        headers = {"User-Agent": "KH-devel"}
        req = urllib.request.Request(url, headers=headers)

        with urllib.request.urlopen(req) as response:
            res = json.loads(response.read())

        # convert to a dict, indexed by pmid
        for x in res:
            pmid = int(x["id"].replace("pmid:", ""))
            citations[pmid] = x["nlm"]["format"]

    return citations


def cluster_articles(pmid_symbol_comat):
    """Clusters articles by their similarity in genes mentioned"""
    print(f"Clustering {pmid_symbol_comat.shape[0]} articles...")

    # mean center
    pmid_symbol_comat = pmid_symbol_comat.apply(lambda x: x - x.mean())

    # apply sign fun and cluster
    pmid_symbol_comat = np.sign(pmid_symbol_comat)

    # choose number of clusters to use for k-means (min: 2)
    # TODO: add parameters as options to query API..
    num_clusters = max(2, int(pmid_symbol_comat.shape[0] / 20))

    kmeans = KMeans(n_clusters=num_clusters, random_state=0).fit(pmid_symbol_comat)
    clusters = kmeans.labels_.tolist()

    res = {}

    for i, pmid in enumerate(pmid_symbol_comat.index):
        res[pmid] = clusters[i]

    return res


def build_network(pmid_symbol_comat):
    """constructs a simple network based on the similarity of article genes mentions"""
    print("Building network..")

    # measure pairwise article correlation
    cor_df = pmid_symbol_comat.T.corr()

    cor_mat = cor_df.to_numpy()
    np.fill_diagonal(cor_mat, 0)

    # for now, keep the top 10% of correlations..
    cutoff = np.quantile(np.abs(cor_mat), 0.9)
    cor_mat[abs(cor_mat) < cutoff] = 0

    ind = np.where(cor_mat > 0)

    pmids = pmid_symbol_comat.index.tolist()

    indFrom = ind[0].tolist()
    indTo = ind[1].tolist()

    edges = []

    for i in range(len(indFrom)):
        article1 = pmids[indFrom[i]]
        article2 = pmids[indTo[i]]

        # number of genes shared by articles
        num_shared = (pmid_symbol_comat.loc[article1] & pmid_symbol_comat.loc[article1]).sum()

        edges.append(
            {
                "source": article1,
                "target": article2,
                "cor": cor_mat[indFrom[i], indTo[i]].item(),
                "num_shared": num_shared.item(),
            }
        )

    return {"links": edges}


@app.post("/query/")
async def query(genes: str, pvalues: Optional[str] = None, disease: Optional[str] = None, 
                max_articles: Optional[int] = 100):
    # split list of genes
    target_symbols = genes.split(",")

    # validate disese MeSH term, if specified
    if disease is not None:
        print(disease)
        if disease not in mesh_pmids.keys():
            return({"error": "Unable to find specified disease MeSH ID"})

    # validate p-value input length
    if pvalues is not None:
        if len(pvalues.split(",")) != len(target_symbols):
            return({"error": "Gene symbol and P-value inputs have different lengths!"})

    # get binary <pmid x gene> matrix
    pmid_symbol_comat = get_pmid_symbol_comat(target_symbols, disease)

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
    # this way, articles with a large number of matched genes are prioritized,
    # penalizing articles which simply include a large number of genes
    #
    article_scores = article_scores** 2 * ratio_matched

    article_scores = article_scores.sort_values(ascending=False).head(max_articles)

    # get subset of article x gene matrix corresponding to the top N hits and cluster
    # articles based on similarity in the genes mentioned..
    pmid_symbol_comat_top = pmid_symbol_comat[pmid_symbol_comat.index.isin(article_scores.index)]

    # drop genes which are not present in any articles, after filtering..
    mask = pmid_symbol_comat_top.sum() > 0

    dropped_genes = ", ".join(sorted(pmid_symbol_comat_top.columns[~mask]))

    print(
        f"Dropping genes that don't appear in any of the top-scoring articles: {dropped_genes}"
    )
    pmid_symbol_comat_top = pmid_symbol_comat_top.loc[:, mask]

    clusters = cluster_articles(pmid_symbol_comat_top)

    # generate network
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

    print("Finalizing results...")

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

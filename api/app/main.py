#!/bin/env python
"""
Pubtator x Gene Search
KH (Oct2021)

- [ ] cluster summary (genes assoc. with each cluster?)
- [ ] cache pmid citation lookups to avoid re-querying server?..
- [ ] add k-means "k" param to api?

"""
import datetime
import json
import numpy as np
import pandas as pd
import urllib.request
from collections import OrderedDict
from biothings_client import get_client
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from sklearn.cluster import KMeans

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

# gene -> pubmed id dataset (entrez/pmids)
gene_pmid_infile = "/data/pmids/gene-pmids.json"

with open(gene_pmid_infile) as fp:
    gene_pmids = json.loads(fp.read())

# pmid -> # unique gene mentions
gene_counts_infile = "/data/gene-counts/pmid_gene_counts.feather"

gene_counts = pd.read_feather(gene_counts_infile).set_index("pmid")


def get_gene_pmid_mat(symbols):
    """Generates an <article x gene> matrix for genes of interest"""
    # convert gene symbols to entrez ids used in pubtator dataset
    mg = get_client("gene")

    mapping = mg.querymany(
        symbols, scopes="symbol", species=9606, fields="entrezgene", as_dataframe=True
    )
    mapping = mapping[~mapping.entrezgene.isna()]

    # check for any genes that could not be mapped and remove them..
    symbols = pd.Series(symbols)

    mask = symbols.isin(mapping.index)

    # exclude genes that could not be mapped, or whose mapped entrez ids are
    # not present in the pubtator dataset
    if mask.sum() < symbols.shape[0]:
        symbols = symbols[mask]

    missing_symbols = mapping.index[~mapping.entrezgene.isin(gene_pmids)]

    if len(missing_symbols) > 0:
        symbols = symbols[~symbols.isin(missing_symbols)]

    target_entrez = mapping.loc[symbols].entrezgene.values

    # subset pmids to get only genes of interest
    dat = {key: gene_pmids[key] for key in target_entrez}

    # get a list of all pubmed ids associated with at least one of the genes
    all_pmids = sorted(set(sum([dat[key] for key in dat], [])))

    print(
        f"Found {len(all_pmids)} articles with one or more of the genes of interest.."
    )

    # convert to a <pmid x gene> matrix
    res: dict[str, list] = {}

    for i, gene in enumerate(target_entrez):
        res[gene] = []

        print(f"Adding entrez gene {gene} ({i + 1}/{len(target_entrez)})...")

        for pmid in all_pmids:
            res[gene].append(pmid in dat[gene])

    dat = pd.DataFrame.from_dict(res)
    dat.columns = symbols
    dat.index = all_pmids

    return dat


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


def cluster_articles(dat):
    """Clusters articles by their similarity in genes mentioned"""
    print(f"Clustering {dat.shape[0]} articles...")

    # mean center
    dat = dat.apply(lambda x: x - x.mean())

    # apply sign fun and cluster
    dat = np.sign(dat)

    # choose number of clusters to use for k-means (min: 2)
    # TODO: add parameters as options to query API..
    num_clusters = max(2, int(dat.shape[0] / 20))

    kmeans = KMeans(n_clusters=num_clusters, random_state=0).fit(dat)
    clusters = kmeans.labels_.tolist()

    res = {}

    for i, pmid in enumerate(dat.index):
        res[pmid] = clusters[i]

    return res


def build_network(dat):
    """constructs a simple network based on the similarity of article genes mentions"""
    print("Building network..")

    # measure pairwise article correlation
    cor_df = dat.T.corr()

    cor_mat = cor_df.to_numpy()
    np.fill_diagonal(cor_mat, 0)

    # for now, keep the top 10% of correlations..
    cutoff = np.quantile(np.abs(cor_mat), 0.9)
    cor_mat[abs(cor_mat) < cutoff] = 0

    ind = np.where(cor_mat > 0)

    pmids = dat.index.tolist()

    indFrom = ind[0].tolist()
    indTo = ind[1].tolist()

    edges = []

    for i in range(len(indFrom)):
        article1 = pmids[indFrom[i]]
        article2 = pmids[indTo[i]]

        # number of genes shared by articles
        num_shared = (dat.loc[article1] & dat.loc[article1]).sum()

        edges.append(
            {
                "source": article1,
                "target": article2,
                "cor": cor_mat[indFrom[i], indTo[i]].item(),
                "num_shared": num_shared.item(),
            }
        )

    return {"links": edges}


@app.get("/query/")
async def query(genes: str, limit: int = 100):
    # split list of genes
    target_symbols = genes.split(",")

    # get <pmid x gene> matrix
    dat = get_gene_pmid_mat(target_symbols)

    # devel: store pmid x gene mat..
    now = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")
    dat.reset_index().to_feather(f"/data/gene-search/{now}.feather")

    # get number of matched genes for each article
    dat_sums = dat.sum(axis=1)

    # score = <matched>^2 * <matched/total>
    # this way, articles with a large number of matched genes are prioritized,
    # penalizing articles which simply include a large number of genes
    pmid_gene_counts = gene_counts.loc[dat_sums.index].num

    scores = dat_sums ** 2 * (dat_sums / pmid_gene_counts)
    scores = scores.sort_values(ascending=False).head(limit)

    # get subset of article x gene matrix corresponding to the top N hits and cluster
    # articles based on similarity in the genes mentioned..
    dat_subset = dat[dat.index.isin(scores.index)]

    # convert to binary
    dat_subset.replace({False: 0, True: 1}, inplace=True)

    # drop genes which are not present in any articles, after filtering..
    mask = dat_subset.sum() > 0

    dropped_genes = ", ".join(sorted(dat_subset.columns[~mask]))

    print(
        f"Dropping genes that don't appear in any of the top-scoring articles: {dropped_genes}"
    )
    dat_subset = dat_subset.loc[:, mask]

    clusters = cluster_articles(dat_subset)

    # generate network
    net = build_network(dat_subset)

    # query pubmed api for article info
    citations = query_article_info(list(dat_subset.index))

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
        clust_submat = dat_subset[dat_subset.index.isin(clust_pmids)]

        clust_counts = clust_submat.sum()
        count_dict = clust_counts.sort_values(ascending=False).to_dict(into=OrderedDict)

        clust_info[clust] = {
            "color": cmap[clust],
            "genes": count_dict,
            "num_articles": len(clust_pmids),
        }

    scores_dict = scores.to_dict()

    dat_subset_sums = dat_sums.loc[dat_subset.index]

    # generate result json, including gene ids, etc. for each article
    articles = []

    for pmid, num_genes in dat_subset_sums.items():

        # retrieve the associated gene names
        pmid_genes = list(dat.loc[pmid][dat.loc[pmid]].index)

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

    res = {"clusters": clust_info, "network": net, "info": info}

    # testing..
    with open(f"/data/gene-search/{now}.json", "w") as fp:
        json.dump(res, fp)

    print("Done!")

    return res

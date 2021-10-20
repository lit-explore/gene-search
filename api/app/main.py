#!/bin/env python
"""
Pubtator x Gene Search
KH (Oct2021)

- [ ] cluster summary (genes assoc. with each cluster?)
- [ ] cache pmid citation lookups to avoid re-querying server?..
- [ ] add k-means "k" param to api?

"""
import json
import numpy as np
import pandas as pd
import urllib.request
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

#  fp = open(gene_pmid_infile)
#  gene_pmids = json.loads(fp.read())
#  fp.close()


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

    if mask.sum() < symbols.shape[0]:
        # missing = symbols[~mask].tolist()
        symbols = symbols[mask]

    target_entrez = mapping.loc[symbols].entrezgene.values

    # exclude any genes with no associated entrez id in the gene/pmid data
    target_entrez = [x for x in target_entrez if x in gene_pmids]

    # subset pmids to get only genes of interest
    dat = {key: gene_pmids[key] for key in target_entrez}

    # get a list of all pubmed ids associated with at least one of the genes
    all_pmids = sorted(set(sum([dat[key] for key in dat], [])))

    print(
        f"Found {len(all_pmids)} articles with one or more of the genes of interest.."
    )

    # convert to a <pmid x gene> matrix
    res: dict[str, list] = {}

    for gene in target_entrez:
        res[gene] = []

        print(gene)

        for pmid in all_pmids:
            res[gene].append(pmid in dat[gene])

    dat = pd.DataFrame.from_dict(res)
    dat.columns = symbols
    dat.index = all_pmids

    return dat


def query_article_info(pmids):
    """Uses PubMed API to retrieve additional information for each article hit"""
    base_url = "https://api.ncbi.nlm.nih.gov/lit/ctxp/v1/pubmed/?format=citation&id="

    url = base_url + ",".join([str(x) for x in pmids])

    headers = {"User-Agent": "KH-devel"}
    req = urllib.request.Request(url, headers=headers)

    with urllib.request.urlopen(req) as response:
        res = json.loads(response.read())

    # convert to a dict, indexed by pmid
    citations = {int(x["id"].replace("pmid:", "")): x["nlm"]["format"] for x in res}

    return citations


def cluster_articles(dat):
    """Clusters articles by their similarity in genes mentioned"""
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
    cor_df = dat.T.corr()

    cor_mat = cor_df.to_numpy()
    np.fill_diagonal(cor_mat, 0)

    # for now, keep the top 4% of correlations..
    cutoff = np.quantile(np.abs(cor_mat), 0.96)
    cor_mat[abs(cor_mat) < cutoff] = 0

    ind = np.where(cor_mat > 0)

    pmids = dat.index.tolist()

    indFrom = ind[0].tolist()
    indTo = ind[1].tolist()

    edges = []

    for i in range(len(indFrom)):
        edges.append((pmids[indFrom[i]], pmids[indTo[i]]))

    return edges


@app.get("/query/")
async def query(genes: str, limit: int = 100):
    # split list of genes
    target_symbols = genes.split(",")

    # get <pmid x gene> matrix
    dat = get_gene_pmid_mat(target_symbols)

    # get gene counts for each article, and sort from highest to lowest
    dat_summary = dat.sum(axis=1).sort_values(ascending=False).head(limit)

    # get subset of article x gene matrix corresponding to the top N hits and cluster
    # articles based on similarity in the genes mentioned..
    dat_subset = dat[dat.index.isin(dat_summary.index)]

    # convert to binary
    dat_subset.replace({False: 0, True: 1}, inplace=True)

    # exclude genes mentioned in very few articles in the subset..
    min_cutoff = dat_subset.shape[0] / 20
    dat_subset = dat_subset.loc[:, dat_subset.sum() >= min_cutoff]

    clusters = cluster_articles(dat_subset)

    # generate network
    net = build_network(dat_subset)

    # query pubmed api for article info
    citations = query_article_info(list(dat_summary.index))

    # generate result json, including gene ids, etc. for each article
    articles = []

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

    for pmid, num_genes in dat_summary.items():

        # retrieve the associated gene names
        pmid_genes = list(dat.loc[pmid][dat.loc[pmid]].index)

        articles.append(
            {
                "id": pmid,
                "genes": ", ".join(pmid_genes),
                "num_genes": num_genes,
                "citation": citations[pmid],
                "cluster": clusters[pmid],
                "color": cmap[clusters[pmid]],
            }
        )

    res = {"articles": articles, "network": net}

    return res

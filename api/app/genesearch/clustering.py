import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

def cluster_articles(pmid_symbol_comat_df:pd.DataFrame):
    """Clusters articles by their similarity in genes mentioned"""
    print(f"Clustering {pmid_symbol_comat_df.shape[0]} articles...")

    # mean center
    pmid_symbol_comat = pmid_symbol_comat_df.apply(lambda x: x - x.mean())

    # apply sign function and cluster
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



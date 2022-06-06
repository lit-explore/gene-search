import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

def cluster_articles(comat:pd.DataFrame) -> dict:
    """
    Clusters articles by their similarity in genes mentioned

    Arguments
    ---------
    comat: pd.DataFrame
        article x gene co-occurrence matrix

    Return
   -------
   out: dictionary mapping from article pmids -> cluster ids
    """
    # mean center and convert to binary
    comat = np.sign(comat.apply(lambda x: x - x.mean()))

    # choose number of clusters to use for k-means (min: 2)
    num_clusters = max(2, int(comat.shape[0] / 20))

    kmeans = KMeans(n_clusters=num_clusters, random_state=0).fit(comat)
    clusters = kmeans.labels_.tolist()

    res = {}

    for i, pmid in enumerate(comat.index):
        res[pmid] = clusters[i]

    return res

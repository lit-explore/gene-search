import numpy as np
import pandas as pd

def build_network(comat:pd.DataFrame) -> dict:
    """
    Constructs a simple network based on the similarity of article genes mentions

    Arguments
    ---------
    comat: pd.DataFrame
        article x gene co-occurrence matrix
    """
    # measure pairwise article correlation
    cor_df = comat.T.corr()

    cor_mat = cor_df.to_numpy()
    np.fill_diagonal(cor_mat, 0)

    # for now, keep the top 10% of correlations..
    cutoff = np.quantile(np.abs(cor_mat), 0.9)
    cor_mat[abs(cor_mat) < cutoff] = 0

    ind = np.where(cor_mat > 0)

    pmids = comat.index.tolist()

    indFrom = ind[0].tolist()
    indTo = ind[1].tolist()

    edges = []

    for i in range(len(indFrom)):
        article1 = pmids[indFrom[i]]
        article2 = pmids[indTo[i]]

        # number of genes shared by articles
        num_shared = (comat.loc[article1] & comat.loc[article1]).sum()

        edges.append(
            {
                "source": article1,
                "target": article2,
                "cor": cor_mat[indFrom[i], indTo[i]].item(),
                "num_shared": num_shared.item(),
            }
        )

    return {"links": edges}



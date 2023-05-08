"""
Pubmed-related functions
"""
import logging
from typing import List, Optional
import numpy as np
import pandas as pd
import json
import sys
import urllib.request

def build_article_gene_comat(entrez_ids:List[int],
                             article_mapping:dict,
                             pmids_allowed: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Generates an <article x gene> matrix for genes of interest

    Arguments
    ---------
    entrez_ids: list[int]
        List of target entrez gene ids
    article_mapping: dict
        Dictionary mapping from entrez gene ids to pubmed article ids
    pmids_allowed: list[str]
        Optional list of PMIDs to restrict result to

    Return
    ------
    out: pd.DataFrame
        Gene x Pubmed article co-occurrence matrix as a pandas dataframe
    """
    # setup logger
    logging.basicConfig(stream=sys.stdout,
                        format='%(asctime)s [%(levelname)s] %(message)s',
                        level=logging.INFO)
    logger = logging.getLogger('gene-search')

    # get a list of all unique pubmed article ids associated with >= 1 gene of interest
    # TODO (May 21, 2022): add optional step here to limit result to articles associated
    # with >= N genes?
    pmid_lists = [article_mapping[entrez_id] for entrez_id in article_mapping]
    matching_pmids = pd.Series(sorted(set(sum(pmid_lists, []))))

    num_matches = len(matching_pmids)
    logger.info("Found %d articles with one or more of the genes of interest..", num_matches)

    # if disease specified, filter to include only articles relating to that disease
    if pmids_allowed is not None:
        num_before = num_matches
        matching_pmids = matching_pmids[matching_pmids.isin(pmids_allowed)]
        num_after = len(matching_pmids)

        logger.info("Excluding %d/%d articles not included in allowed PMIDs.",
                    num_before - num_after, num_before)

    # iterate over genes and create binary vectors corresponding to
    # their presence/absence in each article
    rows = []

    for _, gene in enumerate(entrez_ids):
        row = matching_pmids.isin(article_mapping[gene]).astype(np.int64)
        rows.append(row)

    article_gene_mat = pd.concat(rows, axis=1)
    article_gene_mat.columns = pd.Index(entrez_ids)
    article_gene_mat.index = pd.Index(matching_pmids)

    return article_gene_mat

def query_article_info(pmids:list[str]) -> dict:
    """
    Uses PubMed API to retrieve additional information for each article hit

    Arguments
    ---------
    pmids: list[str]
        List of pubmed article ids
    """

    # setup logger
    logging.basicConfig(stream=sys.stdout,
                        format='%(asctime) [%(levelname)s] %(message)s',
                        level=logging.INFO)
    logger = logging.getLogger('gene-search')

    base_url = "https://api.ncbi.nlm.nih.gov/lit/ctxp/v1/pubmed/?format=citation&id="

    citations = {}

    # helper function chunk a list
    # https://stackoverflow.com/a/312464/554531
    def chunks(lst, n):
        for i in range(0, len(lst), n):
            yield lst[i : i + n]

    pmid_chunks = list(chunks(pmids, 50))

    for i, chunk in enumerate(pmid_chunks):
        logger.info("Querying PubMed API for citation info.. (%d/%d)", i + 1, len(pmid_chunks))

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

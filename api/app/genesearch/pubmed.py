"""
Pubmed-related functions
"""
from typing import List, Optional
import numpy as np
import pandas as pd
import json
import urllib.request

def get_pmid_symbol_comat(target_symbols:List[str], gene_mapping:pd.DataFrame, 
                          entrez_pmids:dict, pmids_allowed: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Generates an <article x gene> matrix for genes of interest

    Arguments
    ---------
    target_symbols: list[str]
        List of gene target gene symbols
    pmids_allowed: list[str]
        Optional list of PMIDs to restrict result to

    Return
    ------
    out: pd.DataFrame
        Gene x Pubmed article co-occurrence matrix as a pandas dataframe
    """
    # check for any genes that could not be mapped and remove them..
    target_series = pd.Series(target_symbols)

    # exclude genes that could not be mapped, or whose mapped entrez ids are
    # not present in the pubtator dataset
    mask = target_series.isin(gene_mapping.index)
    target_series = target_series[mask]

    # convert gene target_series to entrez ids (entrez_pmids is indexed by string ids)
    target_entrez = gene_mapping.loc[target_series].entrezgene.values
    target_entrez = [str(x) for x in target_entrez]

    # exlude genes for which no entrez is present in the entrez -> pmid mapping
    mask = pd.Series([x in entrez_pmids.keys() for x in target_entrez])
    target_entrez = pd.Series(target_entrez)[mask].values

    target_series = target_series[mask.values]

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
    if pmids_allowed is not None:
        num_before = num_matches
        matching_pmids = matching_pmids[matching_pmids.isin(pmids_allowed)]
        num_after = len(matching_pmids)

        print(f"Excluding {num_before - num_after}/{num_before} articles not included in allowed PMIDs.")

    # convert to a binary <pmid x gene> co-occurrence matrix
    pmid_symbol_rows = []

    for i, gene in enumerate(target_entrez):
        row = matching_pmids.isin(entrez_pmids_subset[gene]).astype(np.int64)
        pmid_symbol_rows.append(row)

    pmid_symbol_mat = pd.concat(pmid_symbol_rows, axis=1)

    # convert boolean to numeric (0|1)
    #  pmid_symbol_mat.replace({False: 0, True: 1}, inplace=True)

    pmid_symbol_mat.columns = target_series
    pmid_symbol_mat.index = matching_pmids

    return pmid_symbol_mat

def query_article_info(pmids:list[str]):
    """
    Uses PubMed API to retrieve additional information for each article hit

    Arguments
    ---------
    pmids: list[str]
        List of pubmed article ids
    """

    base_url = "https://api.ncbi.nlm.nih.gov/lit/ctxp/v1/pubmed/?format=citation&id="

    citations = {}

    # helper function chunk a list
    # https://stackoverflow.com/a/312464/554531
    def chunks(lst, n):
        for i in range(0, len(lst), n):
            yield lst[i : i + n]

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

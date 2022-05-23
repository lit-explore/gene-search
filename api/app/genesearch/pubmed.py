"""
Pubmed-related functions
"""
from typing import List, Optional
import numpy as np
import pandas as pd
import json
import urllib.request

def get_pmid_symbol_comat(target_symbols:List[str], gene_mapping:pd.DataFrame, 
                          pubtator_entrez_pmids:dict, pmids_allowed: Optional[List[str]] = None) -> pd.DataFrame:
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
    target_genes = pd.Series(target_symbols)

    # exclude any genes that either can't be mapped to entrez ids
    mask = target_genes.isin(gene_mapping.index)
    target_genes = target_genes[mask]

    # convert gene symbols to entrez ids
    target_entrez = gene_mapping.loc[target_genes].entrezgene.values
    target_entrez = pd.Series([str(x) for x in target_entrez])

    # exclude any genes that aren't present in the pubtator data
    #mask = pd.Series([x in pubtator_entrez_pmids.keys() for x in target_entrez])
    mask = target_entrez.isin(pubtator_entrez_pmids.keys())

    #  target_entrez = pd.Series(target_entrez)[mask].values
    valid_entrez_ids = target_entrez[mask].values

    target_genes = target_genes[mask.values]

    # get subset of gene/article mapping for genes of interest
    entrez_pmids_subset = {}

    for entrez_id in valid_entrez_ids:
        entrez_pmids_subset[entrez_id] = pubtator_entrez_pmids[entrez_id]
        entrez_pmids_subset[entrez_id]

    # get a list of all unique pubmed article ids associated with >= 1 gene of interest
    # TODO (May 21, 2022): add optional step here to limit result to articles associated
    # with >= N genes?
    pmid_lists = [entrez_pmids_subset[entrez_id] for entrez_id in entrez_pmids_subset]
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

    for i, gene in enumerate(valid_entrez_ids):
        row = matching_pmids.isin(entrez_pmids_subset[gene]).astype(np.int64)
        pmid_symbol_rows.append(row)

    pmid_symbol_mat = pd.concat(pmid_symbol_rows, axis=1)

    # convert boolean to numeric (0|1)
    #  pmid_symbol_mat.replace({False: 0, True: 1}, inplace=True)

    pmid_symbol_mat.columns = target_genes
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

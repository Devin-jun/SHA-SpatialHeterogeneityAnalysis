#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import scanpy as sc
import read 
import anndata 
import spatialdata



def cluster_cells(
    s:  spatialdata,
    n_neighbors: int,
    n_pcs: int,
) -> spatialdata:
    """
    clustering cells
    
    Parameters
    ----------
    s
        a spatialdata class
    n_neighbors
        n x neighbors
    n_pcs
        n x pcs
    """
    x = s.E
    df = pd.DataFrame(x.values.T, index=x.columns, columns=x.index)
    adata = anndata.AnnData(X=df)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca_variance_ratio(adata, log=True)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.leiden(adata)
    obs = adata.obs
    s.obs = pd.DataFrame(obs.values.T, index=obs.columns,  columns=obs.index)
    return s



def find_marker_genes(
    s: spatialdata,
    n_genes: int,
) -> dict:
    """
    Finding marker genes by wilcoxon
    
    Parameters
    ----------
    s
        a spatialdata class
    n_genes
        n x top highly differential genes
    """
    x = s.E
    df = pd.DataFrame(x.values.T, index=x.columns, columns=x.index)
    obs = pd.DataFrame(s.obs.values.T, index=s.obs.columns,  columns=s.obs.index)
    adata = anndata.AnnData(X=df, obs=obs)
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups(adata, n_genes=n_genes, sharey=False)
    return adata.uns['rank_genes_groups']




def rename_cluster(
    s: spatialdata,
    newnames: list,
) -> spatialdata:
    """
    mark the cell types
    
    Parameters
    ----------
    s
        a spatialdata class
    newnames:
        newname of each cluster
    """
    types = []
    for index in s.obs.index:
        for n in range(len(s.obs.loc[index].values)):
            cluster = int(s.obs.loc[index].values[n])
            s.obs.loc[index].values[n] = newnames[cluster]
    return s


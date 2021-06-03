import pandas as pd
from sklearn.metrics.pairwise import pairwise_distances
import math


def weight_matrix(cor, ann):
    """
    construct weight matrix for laplacian score calculation
    Parameters
    ----------
    cor：cor matrix of cells
    ann: annotations of types of cells

    Returns
    -------
    normalization weight matrix
    """
    n_length = len(cor)
    d_matrix = pairwise_distances(cor)
    d_max = d_matrix.max()
    d_min = d_matrix.min()
    for i in range(n_length):
        for j in range(n_length):
            if i == j:
                d_matrix[i][j] = 0
            elif ann[i] == ann[j]:
                d_matrix[i][j] = math.exp((d_matrix[i][j] - d_min) / (d_max - d_min))
            else:
                d_matrix[i][j] = 0
    return d_matrix


def laplacian_score(w_matrix, e_matrix, sort='decending'):
    """
    calculate laplacian score for every
    Parameters
    ----------
    w_matrix：weight matrix of cells
    e_matrix: expression matrix of single cell
    sort：ascending or descending

    Returns
    -------
    laplacian score dataframe
    """
    n_genes = e_matrix.shape[0]
    n_cells = e_matrix.shape[1]
    feature_dict = {}
    for i in range(n_genes):
        name = e_matrix.iloc[i].name
        var_feature = e_matrix.iloc[i].var()
        fr = 0.0
        for cell1 in range(n_cells):
            for cell2 in range(cell1 + 1, n_cells):
                fr += ((e_matrix.iloc[i][cell1] - e_matrix.iloc[i][cell2]) ** 2) * w_matrix[cell1][cell2]
        lp_score = fr / var_feature
        feature_dict[name] = lp_score
        feature_df = pd.DataFrame.from_dict(feature_dict, orient='index', columns=['lp_score'])
        feature_df.index.name = 'genes'
        if sort == 'decending':
            sorted_df = feature_df.sort_values(by='lp_score', ascending=False)
        else:
            sorted_df = feature_df.sort_values(by='lp_score', ascending=True)
    return sorted_df

def filtergenes(n_genes,exp):
    """

    Parameters
    ----------
    n_genes: the Nth genes with biggest variance
    exp: expression matrix

    Returns
    -------
    the names of Nth genes
    """
    genes_len = len(exp)
    genes_var = {}
    for i in range(genes_len):
        genes_var[exp.iloc[i].name] = exp.iloc[i].var()
    genes_order = sorted(genes_var.items(), key=lambda x: x[1], reverse=True)
    filter_indexs = [genes_order[j][0] for j in range(n_genes)]
    return pd.DataFrame(exp, index=filter_indexs)


def cell_laplacian_score(celltype, coor, ann, exp, sort='decending'):
    """

    Parameters
    ----------
    celltype：the cell types or the cluster id
    coor: coord matrix of cells
    ann: the cell_types or cluster types of cell
    exp: expression matrix
    sort: decending

    Returns
    -------
    the spatial genes and score
    """

    indexs = [ann.iloc[i].name for i in range(len(ann)) if ann.iloc[i]['cell_types'] == celltype]
    new_coor = pd.DataFrame(coor, index=indexs)
    new_ann = pd.DataFrame(ann, index=indexs)
    new_exp = pd.DataFrame(exp, columns=indexs)
    n_length = len(new_coor)
    d_matrix = pairwise_distances(new_coor)
    d_max = d_matrix.max()
    d_min = d_matrix.min()
    for i in range(n_length):
        for j in range(n_length):
            if i == j:
                d_matrix[i][j] = 0
            else:
                d_matrix[i][j] = math.exp(-(d_matrix[i][j] - d_min) / (d_max - d_min))
    ls = laplacian_score(d_matrix, new_exp, sort)
    return ls
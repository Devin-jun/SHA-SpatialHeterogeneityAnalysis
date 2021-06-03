import numpy as np
import pandas as pd
import read
import spatialdata
import math

#filter cell
def filter_cell(exp,min_genes=None,max_genes=None):
    """
    filter low quality cells
    """
    if min_genes == None or type(min_genes)!=int:
        raise ValueError("min_genes should be int ")
    x_min = exp[exp > 0].count() > min_genes
    exp = exp.iloc[:,np.array(x_min).astype(bool)]
    if type(max_genes)==int:
        x_max = exp[exp > 0].count() < max_genes
        exp = exp.iloc[:,np.array(x_max).astype(bool)]
    return exp
#filter genes

def filter_gene(exp,min_count=None,max_count=None):
    """
    filter low quality genes
    """
    if min_count == None or type(min_count)!=int:
        raise ValueError("min_count should be int ")
    x_min = exp[exp > 0].count(axis=1) > min_count
    exp = exp.iloc[np.array(x_min).astype(bool),:]
    if type(max_count)==int:
        x_max = exp[exp > 0].count(axis=1) < max_count
        exp = exp.iloc[np.array(x_max).astype(bool),:]
    return exp


def sha_norm(exp, c=10e5):
    """
    normalization exp data
    """
    exp_sum = exp.sum(axis=0)
    norm_exp = exp * c / exp_sum
    norm_exp = norm_exp.apply(lambda x: np.log(x + 1))
    norm_exp = norm_exp.apply(lambda x: round(x, 2))

    return norm_exp
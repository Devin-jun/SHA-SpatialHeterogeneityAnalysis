#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import anndata
import csv

def read_csv(
    filename: str,
    delimiter: str = ",",
) -> pd.DataFrame:
    """
    Read `.csv` file.
    Parameters
    ----------
    filename
        Data file.
    delimiter
        Delimiter that separates data within text file.
        If `None`, will split at arbitrary number of white spaces,
        which is different from enforcing splitting at single white space `' '`.
    first_column_names
        Assume the first column stores row names.
    """
    temp = []
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter = delimiter)
        for row in reader:
            temp.append(row)
    df = pd.DataFrame(temp[1:],columns = temp[0])
    return df





def read_mat(
    filename: str,
) -> pd.DataFrame:
    """
    Read `.mtx` file.

    Parameters
    ----------
    filename
        The filename.
    """
    from scipy.io import mmread

    # could be rewritten accounting for dtype to be more performant
    X = mmread(fspath(filename)).astype(dtype)
    from scipy.sparse import csr_matrix
    X = csr_matrix(X)
    df = pd.DataFrame(X)
    return df




def read_txt(
    filename: str,
    delimiter: str = None,
    column_names: list = None,
) -> pd.DataFrame:
    """\
    Read `.txt`(text) file.

    Parameters
    ----------
    filename
        Data file, filename or stream.
    delimiter
        Delimiter that separates data within text file. If `None`, will split at
        arbitrary number of white spaces, which is different from enforcing
        splitting at single white space `' '`.
    first_column_names
        Assume the first column stores row names.
    """
    if column_names:
        df = pd.read_csv(filename, header=0, sep=delimiter, index_col=column_names)
        return df
    df = pd.read_csv(filename, header=0, sep=delimiter)
    return df

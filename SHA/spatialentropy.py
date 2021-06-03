from math import log
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict
from itertools import combinations_with_replacement, product
from typing import Sequence
from sklearn.metrics import pairwise_distances
import spatialdata


def reduce_matrix_size(matrix, types, dicts):
    """
    reduce matrix size
    Parameters
    ----------
    matrix: array
    types: array
    storage: dict

    Returns
    -------
    array,row_merged matrix based on types
    """
    for i, arr in enumerate(matrix):
        dicts[types[i]].append(arr)

    for k, v in dicts.items():
        dicts[k] = np.asarray(v).sum(axis=0)

    new_types = []
    new_matx = []
    for k, v in dicts.items():
        new_types.append(k)
        new_matx.append(v)

    return new_matx

def type_adj_matrix(matrix, types):
    """return an N * N matrix, N is the number of unique types

    Args:
        matrix: array
        types: array

    Returns:
         tuple, matrix and the unique types

    """
    unitypes = np.unique(types)

    storage = OrderedDict(zip(unitypes, [[] for _ in range(len(unitypes))]))
    new_matrix = reduce_matrix_size(matrix, types, storage)

    storage = OrderedDict(zip(unitypes, [[] for _ in range(len(unitypes))]))
    type_matrix = reduce_matrix_size(np.asarray(new_matrix).T, types, storage)
    return np.array(type_matrix), unitypes

def pairs_counter(matrix, types, order=False):
    """count how many pairs of types in the matrix

    Args:
    matrix: array
    types: array
    order: bool, if True, (x1, x2) and (x2, x1) is not the same

    Returns:
        dict, the count of each pairs

    """
    it = np.nditer(matrix, flags=["multi_index"])

    if order:
        combs = [i for i in product(types, repeat=2)]
        storage = OrderedDict(zip(combs, [0 for _ in range(len(combs))]))

        for x in it:
            (i1, i2) = it.multi_index
            storage[(types[i1], types[i2])] += x
    else:
        combs = [i for i in combinations_with_replacement(types, 2)]
        storage = OrderedDict(zip(combs, [0 for _ in range(len(combs))]))

        for x in it:
            (i1, i2) = it.multi_index
            if i1 <= i2:
                storage[(types[i1], types[i2])] += x
            else:
                storage[(types[i2], types[i1])] += x

    return storage

def interval_pairs(arr):
    new_arr = []
    for i, x in enumerate(arr):
        if i < len(arr) - 1:
            new_arr.append((x, arr[i + 1]))

    return new_arr

def cell_ShannonEntropy(
        obs: pd.DataFrame,
):
    """

    Parameters
    ----------
    obs
        Annotation  data matrix of M cells
    -------
    calculate the spatial shannon entropy of each cell type
    """
    n_cell_types = set(obs['cell_types'])
    cell_shannon = {}
    for i in n_cell_types:
        cell_type = obs.iloc[list(np.where(obs['cell_types'] == i)[0])]
        length_cell_type = len(cell_type)
        cell_type_dict = dict(cell_type['label'].value_counts())
        shannon_ent = 0.0
        for key in cell_type_dict:
            prob = float(cell_type_dict[key]) / length_cell_type
            shannon_ent -= prob * log(prob, 2)
        cell_shannon[i] = shannon_ent

    return cell_shannon


class spatial_entropy(object):
    """
    calculate the spatial entropy  by altieri entropy
    """

    def __init__(self, cell_cor, types, cut=None, order=False, base=None):
        if len(cell_cor) != len(types):
            raise ValueError("length of cell and cell types should be the same")
        if base is None:
            base = np.e
        self._cell_cor = cell_cor
        self._types = types
        self._order = order
        self._base = base
        self.adj_matrix = pairwise_distances(self._cell_cor)

        if isinstance(cut, int):
            self._break = interval_pairs(np.linspace(0, self.adj_matrix.max(), cut + 2))
        elif isinstance(cut, Sequence):
            self._break = interval_pairs(cut)

        elif cut is None:
            self._break = interval_pairs(np.linspace(0, self.adj_matrix.max(), 3))  # 没有指定cut长度则最小划分为三段

        else:
            raise ValueError("'cut' must be an int or an array-like object")

        self._wrap()

    def _Z_W(self):

        zw = []
        for (p1, p2) in self._break:
            bool_matx = ((self.adj_matrix > p1) & (self.adj_matrix <= p2)).astype(int)  # bool矩阵计算出在p1和p2区间内的坐标点
            type_matx, utypes = type_adj_matrix(bool_matx, self._types)
            pairs_counts = pairs_counter(type_matx, utypes, self._order)
            zw.append(pairs_counts)

        return zw

    def _Z(self):

        bool_matx = (self.adj_matrix >= 0).astype(int)
        type_matx, utypes = type_adj_matrix(bool_matx, self._types)
        z = pairs_counter(type_matx, utypes, self._order)

        return z

    def _W(self):

        w = []
        for (p1, p2) in self._break:
            w.append(p2 - p1)

        w = np.asarray(w)

        w = w / w.sum()

        return w

    def _wrap(self):

        zw = np.asarray(self._Z_W())
        z = self._Z()
        w = np.asarray(self._W())

        pz = np.array(list(z.values()))
        pz = pz / pz.sum()

        H_Zwk = []  # H(Z|w_k)
        PI_Zwk = []  # PI(Z|w_k)

        for i in zw:
            v_ = i.values()

            v, pz_ = [], []
            for ix, x in enumerate(v_):
                if x != 0:
                    v.append(x)
                    pz_.append(pz[ix])

            v = np.asarray(v)
            pz_ = np.asarray(pz_)

            v = v / v.sum()
            H = v * np.log(1 / v) / np.log(self._base)
            PI = v * np.log(v / pz_) / np.log(self._base)
            H_Zwk.append(H.sum())
            PI_Zwk.append(PI.sum())

        self.residue = (w * np.asarray(H_Zwk)).sum()
        self.mutual_info = (w * np.asarray(PI_Zwk)).sum()
        self.entropy = self.mutual_info + self.residue


# 指定单个区域的空间熵计算
def swspatial_entropy(df, ux, uy, l, w, span, d='h', w_cut=10):
    """

    Parameters
    ----------
    data: the dataframe contain coordinates of each cell and cluster id 
    ux: upleft x of the window,coordinate space
    uy: upleft y of the window,coordiante sapce
    l:length of the window,coordinate space
    w:width of the window,coordinate space
    span: the step size of window
    d: h,horizontal direction or v,vertical direction
    Returns
    -------
    the spatial entropy of the window area
    """""
    if l <= 0 or w <= 0:
        raise ValueError("length and width of the window should be greater than 0")
    if d == 'h':
        x_coord_max = df['x_coord'].max()
        if ux + l <= x_coord_max:
            site = []
            swse = []
            for index in range(ux, x_coord_max, span):
                spot = df.loc[df['x_coord'] > index]
                spot = spot.loc[spot['x_coord'] < (index + l)]
                spot = spot.loc[spot['y_coord'] > uy]
                spot = spot.loc[spot['y_coord'] < uy + w]
                coord = np.array(spot[['x_coord', 'y_coord']])
                clusters = list(spot.iloc[:, 2])
                print(len(coord))
                print(len(clusters))
                se = spatial_entropy(coord, clusters, cut=w_cut)
                swse.append(se)
                site.append(index)
    if d == 'v':
        y_coord_max = df['y_coord'].max()
        if uy + w <= y_coord_max:
            site = []
            swse = []
            for index in range(uy, y_coord_max, span):
                spot = df.loc[df['x_coord'] > ux]
                spot = spot.loc[spot['x_coord'] < ux + l]
                spot = spot.loc[spot['y_coord'] > index]
                spot = spot.loc[spot['y_coord'] < index + w]
                coord = np.array(spot[['x_coord', 'y_coord']])
                clusters = list(spot.iloc[:, 2])
                se = spatial_entropy(coord, clusters, cut=w_cut)
                swse.append(se)
                site.append(index)
    se_df = pd.DataFrame({'ul_of_windows': site, 'spatial_entropy': swse})
    return se_df


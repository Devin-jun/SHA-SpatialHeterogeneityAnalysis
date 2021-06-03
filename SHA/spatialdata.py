#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd



class SpatialData(object):
    """
    An annotated data class especially for single-cell spatial data 
    Parameters 
    -------
    E
         A #genes  × #cells data matrix.
    obs
        Key indexed one-dimensional observations annotation of length #cells. 
    cor 
        Key indexed one-dimensional coordinate (geometry)  of length #cells.
    var 
        Key indexed one-dimensional variables annotation of length #genes.
    dtype
        Data type used for storage.
    shape
        Shape tuple (#genes, #cells). Can only be provided if `E` is `None`.
    """
    
    def __init__(
        self,
        E: pd.DataFrame = None,
        obs: pd.DataFrame = None,
        cor: pd.DataFrame = None,
        var: pd.DataFrame = None,
     ): 
        self.E = E
        self.obs = obs
        self.cor = cor
        self.var = var
        
    @property
    def shape_of_E(self):
        """
        shape of matrix E
        """
        if self.E.empty == False:
            E_row = self.E.shape[0]
            E_col = self.E.shape[1]
            print('E with n_genes x n_cells = {} x {}'.format(E_row, E_col))
            return (E_row, E_col)
        return "E is not available"
    
              


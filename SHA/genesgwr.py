import numpy as np
import pandas as pd
import read
import spatialdata
import math
import pysal.lib as ps
from mgwr.gwr import GWR
import matplotlib.pyplot as plt

def gene_reg_gwr(coords,y,x,bw,fixed,kernel):
    """
    @description  :
    construct the genes regulate netwaork by gwr
    ---------
    @param  :
    coor: the coordinates of each cells 
    y: array n*k, independent variable, exlcuding the constant
    x: array n*k, independent variable, exlcuding the constant
    bw: scalar bandwidth value consisting of either a distance or N
                    nearest neighbors; user specified or obtained using
                    Sel_BW
    fixed: boolean True to include intercept (default) in model and False to exclude
                    intercept.
    kernel: string type of kernel function used to weight observations;
                    available options:
                    'gaussian'
                    'bisquare'
                    'exponential'
    -------
    @Returns  :
    gwr object
    -------
    """
    model = GWR(coords, y, x, bw=bw, fixed=fixed, kernel=kernel)
    result = model.fit()
    return result


def plot_gene_gwr(params, coor_x, coor_y, cm='RdYlBu', t):
    """
    @description  :
    plot the paramters distribution in the cells
    ---------
    @param  : 
    params : the value of the paramter
    coor_x: coor X
    coor_y: coor Y
    t: title of the picture
    -------
    @Returns  :
    -------
    """

    ccm = plt.cm.get_cmap(cm)
    sc = plt.scatter(coor_x, coor_y, c=params, cmap=ccm)
    plt.colorbar(sc)
    plt.title(t)
    plt.show()
    
    

    
    
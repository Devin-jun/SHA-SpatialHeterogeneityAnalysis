import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import pairwise_distances
import math
import matplotlib.pyplot as plt
import seaborn as sns

def dist_weight_matrix(cor,k): 

    """
    @description  : construct the k nearest cells weight matrix with the spatial information
    ---------
    @param  :
                cor: the coordinates of each cells 
                k: the number of nearest cells used to construct the weight matrix 
    -------
    @Returns  :
                a_sum: accessibility score of each cell 
    -------
    """
    
    n_length = len(cor)
    d_matrix = pairwise_distances(cor)
    d_max = d_matrix.max()
    d_min = d_matrix.min()
    for i in range(n_length):
        for j in range(n_length):
            if i == j:
                d_matrix[i][j] = 0
            else:
                d_matrix[i][j] = d_matrix[i][j]
    
    sorted_matrix = np.sort(d_matrix)
    
    sa = {}
    for i in range(n_length):
        a_sum = 0
        for j in range(k+1):
            if sorted_matrix[i][j] == 0:
                a_sum += 0
            else:
                a_sum += (math.exp(1/sorted_matrix[i][j])-1)*100
        sa[cor.iloc[i].name] = a_sum
        
    sa_df = pd.DataFrame.from_dict(sa, orient='index',columns=['a_sum'])
    sa_df.index.name = 'ID'
             
    return sa_df


def plot_spatial_accessibility(a_sum,ann,cor,p_sizes=(20,500),p_alpha=0.5,p_palette='muted',p_height=5):
    """
    @description  : plot the spatial accessibility
    ---------
    @param  : 
                a_sum: accessbility score
                ann: annotation of cell types
                cor: the coordinates of each cells
                p_sizes: size of the point 
                p_alpha: the transplant 
                p_height: the height of the plot

    -------
    @Returns  :
    -------
    """
    a_sum['cell_types'] = ann['cell_types']
    a_sum['X'] = cor['X']
    a_sum['Y'] = cor['Y']
    sns.relplot(x="X", y="Y", hue="cell_types", size="a_sum",height=p_height,sizes=p_sizes,
            alpha=p_alpha, palette=p_palette,
             data=a_sum) 
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import array


def distance_compute(
        cor: pd.DataFrame,
) -> pd.DataFrame:
    """
    computing distance between two cells
    Parameters
    ----------
    cor
        The coordinate (geometry) x,y of cells
    -------
    """
    dis_mat = squareform(pdist(cor, metric='euclidean'))
    return dis_mat


def spatialcluster(
        cor: pd.DataFrame,
        eps: int,
        min_pts: int,
) -> array:
    """
    Spatial clustering by DBSCAN
    Parameters
    ----------
    cor
        distance matrix
    eps
        neighborhood distance
    min_pts
        minPts
    -------
    """
    dis_mat = distance_compute(cor)
    n_cells, dim = cor.shape
    core_index = np.where(np.sum(np.where(dis_mat <= eps, 1, 0), axis=0) >= min_pts)[0]
    cell_labels = np.full((n_cells,), -1)
    cluster_id = 0
    for point in core_index:
        if cell_labels[point] == -1:
            cell_labels[point] = cluster_id
            neighbour_cell = np.where((dis_mat[:, point] <= eps) & (cell_labels == -1))[0]
            seeds = set(neighbour_cell)
            while len(seeds) > 0:
                new_point = seeds.pop()
                cell_labels[new_point] = cluster_id
                new_point_neighbor = np.where(dis_mat[:,new_point] <= eps)[0]
                if len(new_point_neighbor) >= min_pts:
                    for new_point_result in new_point_neighbor:
                        if cell_labels[new_point_result] == -1:
                            seeds.add(new_point_result)
            cluster_id = cluster_id + 1
    return cell_labels

def plot_domain(
        cor: pd.DataFrame,
        labels: array,
        colors=['black', 'blue', 'green', 'yellow', 'red', 'purple', 'orange', 'brown'],
        size=12,
):
    """
    plot the spatial cluster result by DBSCAN
    Parameters
    ----------
    cor
        distance matrix
    labels
        cluster label of every cells
    colors
        colors of clusters
    size
        size  of every point
    -------
    """
    try:
        if len(colors) < len(set(labels)):
           raise ValueError("the number of colors is less than the number of clusters")
    except ValueError as e:
        print("error", repr(e))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(-1,len(set(labels))):
        color = colors[i % len(colors)]
        sub_cluster = cor.iloc[np.where(labels == i)]
        ax.scatter(sub_cluster.iloc[:, 0], sub_cluster.iloc[:, 1], c=color, s=size)
    plt.show()


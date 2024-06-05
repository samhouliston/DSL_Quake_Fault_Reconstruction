import numpy as np
import matplotlib.pyplot as plt
import glasbey
from matplotlib.colors import ListedColormap
from sklearn.decomposition import PCA
import colorcet as cc
import plotly.graph_objs as go


def plot_3Ddataset(X: np.ndarray,
                 labels: np.ndarray,
                 marker_sz: float = 4.,
                 cmap = None):
    
    '''
    Create every 2D combination plot from 3D data, colored according to labels

    Parameters
    -----------
    X: np.ndarray
        The 3D data to plot of shape (n_samples, 3)
        
    labels: np.ndarray
        An array of shape (n_samples,) according to which to points are colored

    marker_sz: float
        The size of the markers

    cmap: Colormap
        Colormap to use for the markers, if None will use a glasbey colormap

    
    Returns
    --------
    fig: matplotlib.figure.Figure
        The figure object
        
    ax: np.ndarray
        An array of shape (2,2) containing the ax objects of the figure

    '''

    if cmap is None:
        cmap = ListedColormap(glasbey.create_palette(palette_size=len(np.unique(labels))))
    
    fig, ax = plt.subplots(2,2)
    fig.set_figwidth(10)
    fig.set_figheight(10)

    ax[0,0].scatter(X[:,0], X[:,1], s=marker_sz, c=labels, cmap = cmap)
    ax[0,1].scatter(X[:,2], X[:,1], s=marker_sz, c=labels, cmap = cmap)
    ax[1,0].scatter(X[:,0], X[:,2]*-1, s=marker_sz, c=labels, cmap = cmap)

    ax[0,0].set_xlabel('x')
    ax[0,0].set_ylabel('y')
    ax[0,1].set_xlabel('z')
    ax[0,1].set_ylabel('y')
    ax[1,0].set_xlabel('x')
    ax[1,0].set_ylabel('z')

    ticks = np.linspace(np.max(X[:,2]), np.min(X[:,2]), 5)
    ticks = ax[1,0].get_yticks()[1:-1]
    ax[1,0].set_yticks(ticks, np.round(-1*ticks, 2))

    ax[1,1].axis('off')

    return fig, ax




def plot_components(X: np.ndarray,
                    labels: np.ndarray,
                    marker_sz: float = 2.,
                    cmap = None,
                    X_pca: np.ndarray = None):
    
    '''
    Plot the 3 principal components of the data as 3 plots

    Parameters
    -----------
    X: np.ndarray
        The data as an array of shape (n_samples, n_features)

    labels: np.ndarray
        An array of labels according to which the data points are colored

    marker_sz: float
        The size of the markers

    cmap: Colormap
        Colormap to use for the markers, default = cet_glasbey_light

    X_pca: np.ndarray
        Data to compute the principal components on, if None use X. Default = None

    
    Returns
    --------
    fig: matplotlib.figure.Figure
        The figure object
        
    ax: np.ndarray
        An array of shape (2,2) containing the ax objects of the figure
    '''

    if cmap is None:
        cmap = ListedColormap(glasbey.create_palette(palette_size=len(np.unique(labels))))

    if X_pca is None:
        X_pca = X

    # project the data onto PCs
    pca = PCA(n_components=X_pca.shape[1], random_state=71)
    pca.fit(X_pca)
    
    X = pca.transform(X)

    # plot PCs
    fig, ax = plot_3Ddataset(X, labels, marker_sz, cmap)

    # rename axes
    ax[0,0].set_xlabel('PC 1')
    ax[0,0].set_ylabel('PC 2')
    ax[0,1].set_xlabel('PC 3')
    ax[0,1].set_ylabel('PC 2')
    ax[1,0].set_xlabel('PC 1')
    ax[1,0].set_ylabel('PC 3')

    return fig, ax



def make_3D_plot(X: np.ndarray, 
                 labels: np.ndarray,
                 title: str = '',
                 marker_sz: float = 3.):
    """
    Create an interactive 3D plot from the data colored according to the labels.

    Parameters
    ----------
    X: np.ndarray
        The data points as an array of shape (n_samples, 3)
    
    labels: np.ndarray
        The labels as an array of shape (n_samples,) according to which to color the data

    title: str
        The title of the figure, default = ''

    marker_sz: int
        The size of the individual markers in the plot, default = 3.


    Returns
    --------
    None
    """

    palette = glasbey.create_palette(palette_size=np.max(np.unique(labels))+1)

    fig = go.Figure()
    
    for i in np.unique(labels):
        # Extract points of the current cluster
        X_curr = X[labels == i]

        # Plot the points
        fig.add_trace(go.Scatter3d(
            x=X_curr[:,0],
            y=X_curr[:,1],
            z=X_curr[:,2],
            mode='markers',
            marker=dict(
                size=marker_sz,
                opacity=0.8,
                color=palette[i]
            ),
            name=f'{i}'
        ))
    
    # add labels
    fig.update_layout(scene=dict(
        xaxis_title='X',
        yaxis_title='Y',
        zaxis_title='Z'),
        title=title
    )

    fig.show()

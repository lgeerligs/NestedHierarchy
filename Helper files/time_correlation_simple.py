import matplotlib.patches as patches
import numpy as np


def time_correlation_simple(ax, corr, GSBSbounds=0):
    # ax: where it is plotted
    # data: 2D matrix, time x voxels
    # GSBS (opt): GSBS object

    # Plot the matrix
    ax.imshow(corr, interpolation='none', vmin=-1,vmax=1)
    ax.set_xlabel('TR')
    ax.set_ylabel('TR')

    bounds = np.where(GSBSbounds > 0)[0]-0.5

    for i in range(len(bounds)-1):
        rect = patches.Rectangle(
            (bounds[i], bounds[i]),
            bounds[i + 1] - bounds[i],
            bounds[i + 1] - bounds[i],
            linewidth=4, edgecolor='w', facecolor='none'
        )
        ax.add_patch(rect)

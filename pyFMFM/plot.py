""" 
Plotting utilities and animation into gif
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from params import *
# from loop import *

def plotheatmap(u_k, k, save):
    # Clear the current plot figure
    plt.clf()
    plt.figure(figsize=(4, 4))

    # Selecting the axis-X making the bottom and top axes False.
    plt.tick_params(axis='x', which='both', bottom=False,
                    top=False, labelbottom=False)
      
    # Selecting the axis-Y making the right and left axes False
    plt.tick_params(axis='y', which='both', right=False,
                    left=False, labelleft=False)
    # Iterating over all the axes in the figure
    # and make the Spines Visibility as False
    for pos in ['right', 'top', 'bottom', 'left']:
        plt.gca().spines[pos].set_visible(False)

    # plt.title(f"State variable at t = {k*dt:.3f}")
    # plt.xlabel("x")
    # plt.ylabel("y")

    # This is to plot u_k (u at time-step k)
    plt.pcolormesh(xv, yv, u_k, cmap=plt.cm.Oranges,shading='gouraud')
    # , vmin=0, vmax=plotmax)
    # plt.colorbar()

    if save:
        filename = 'D/sol.png'
        if os.path.isfile(filename):
           os.remove(filename)

        plt.savefig(filename,bbox_inches='tight',pad_inches=0)
        plt.close()
        return
    else:
        return plt

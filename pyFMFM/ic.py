""" 
This has the IC setup and helper functions
"""

import numpy as np
from params import *
import matplotlib.pyplot as plt

def set_ic(uin, coloring=None):
    uic = np.zeros_like(uin)

    if ic == 'hotend':
        u_top = 100.0
        u_left = 0.
        u_bottom = 0.
        u_right = 0.

        # Set initial edge "boundary" conditions
        uic[:buffer_x, :] = u_left
        uic[(Nx-buffer_x):, :] = u_right

        uic[:, :buffer_y] = u_bottom
        uic[:, (Ny-buffer_y):] = u_top

    elif ic == 'zeros':
        u_top = 0.
        u_left = 0.
        u_bottom = 0.
        u_right = 0.

        # Set initial edge "boundary" conditions
        uic[:buffer_x, :] = u_left
        uic[(Nx-buffer_x):, :] = u_right

        uic[:, :buffer_y] = u_bottom
        uic[:, (Ny-buffer_y):] = u_top

    elif ic == '2dgauss':

        if case == 'inhomogeneous':
            mywidth = Lx/20.
        else:
            mywidth = Lx/4.

        uic = make2DGauss(xv, yv, fwhm = mywidth, center=[Lx/2.,Ly/2.])

        u_top = 0.
        u_left = 0.
        u_bottom = 0.
        u_right = 0.

        # Set initial edge "boundary" conditions
        uic[:buffer_x, :] = u_left
        uic[(Nx-buffer_x):, :] = u_right

        uic[:, :buffer_y] = u_bottom
        uic[:, (Ny-buffer_y):] = u_top

    elif ic == '1dgauss':
        uic = make1DGauss(xv, yv, fwhm = 0.025, center=Lx/2.)

    elif ic == 'florian':
        uic[int(Nx/2.),:] = 1./Ny

    elif ic == 'loopy' and coloring is not None:
        if len(coloring) != Nx:
            #Internal coloring
            j = 0
            for i in range(buffer_x,Nx-buffer_x):
                uic[i,:] = coloring[j]
                j += 1

        elif len(coloring) == Nx:
            #Full coloring
            for i in range(Nx):
                uic[i,:] = coloring[i]

    else:
        raise Exception('No such case')

    return uic

def make2DGauss(xv, yv, fwhm, center):
    """ 2D Gaussian kernel.
    fwhm: full-width-half-maximum, or effective radius.
    """

    x0 = center[0]
    y0 = center[1]

    r2 = (xv-x0)**2 + (yv-y0)**2
    ret = np.exp(-4*np.log(2) * r2 / fwhm**2)

    # plt.clf()
    # plt.pcolormesh(xv,yv,ret, cmap=plt.cm.jet, vmin=0, vmax=plotmax)
    # plt.colorbar()
    # plt.savefig('ic.png')
    # plt.close()

    return ret

def make1DGauss(xv, yv, fwhm, center):
    r2 = (xv-center)**2 
    ret = np.exp(-r2 / fwhm)
    return ret


"""
This has various operations, include the rhs 
and derivative operators (central differencing),
also boundary conditions
"""

import numpy as np
import matplotlib.pyplot as plt
from params import *
from numba import njit 

@njit
def apply_bc(uin):
    usave = uin.copy()

    if bc_x == 'none':
        pass
    elif bc_x == 'periodic':
        uin[0,:] = usave[Nx-2,:]
        uin[Nx-1,:] = usave[1,:]
    elif bc_x == 'noflux':
        uin[0,:] = usave[1,:]
        uin[Nx-1,:] = usave[Nx-2,:]
    else:
        raise Exception('No such x-dir BC, bc_x')

    if bc_y == 'none':
        pass
    elif bc_y == 'periodic':
        uin[:,0] = usave[:,Ny-2]
        uin[:,Ny-1] = usave[:,1]
    elif bc_y == 'noflux':
        uin[:,0] = usave[:,1]
        uin[:,Ny-1] = usave[:,Ny-2]
    else:
        raise Exception('No such y-dir BC, bc_y')

    return uin

@njit
def rhs(uin2):
    myrhs = np.zeros_like(uin2)
    uin = apply_bc(uin2)

    for i in range(buffer_x, Nx-buffer_x):
        for j in range(buffer_y, Ny-buffer_y):
            if (case == 'heat'): 
                myrhs[i,j] = alpha*laplace2d(uin,i,j)
            elif (case == 'advection'): 
                myrhs[i,j] = -adv_speed*dudx(uin,i,j)
            elif (case == 'advection-diffusion'): 
                myrhs[i,j] = 1*dudy(uin,i,j) + du2dy2(uin,i,j)
            elif (case == 'parallel'): 
                myrhs[i,j] = -np.cos(y[j]) * dudx(uin,i,j) \
                    + du2dy2(uin,i,j) 
            elif (case == 'inhomogeneous' or case == 'inhomogeneous-periodic'): 
                myrhs[i,j] = \
                      -1*u1(x[i],y[j]) * dudx(uin,i,j) \
                    + -1*u2(x[i],y[j]) * dudy(uin,i,j) \
                    +  du2dy2(uin,i,j) \
                    + 0.05 * du2dx2(uin,i,j) \
                    + forcing(x[i])
            else:
                raise Exception('No such case')

    return myrhs



@njit
def forcing(xin):
    if case == 'inhomogeneous':
        ret = 1. 
    elif case == 'inhomogeneous-periodic':
        ret = np.cos(xin-np.pi)
    # else:
    #     raise Exception('No such case, u1')

    return ret 


@njit
def u1(xin,yin):
    if case == 'inhomogeneous':
        ret = (1+np.cos(xin-np.pi))*np.cos(yin)
    elif case == 'inhomogeneous-periodic':
        ret = 0.5*(2.+np.cos(xin-np.pi))*np.cos(yin)
    else:
        raise Exception('No such case, u1')

    return ret 

@njit
def u2(xin,yin):
    if case == 'inhomogeneous':
        ret = np.sin(xin-np.pi)*np.sin(yin)
    elif case == 'inhomogeneous-periodic':
        ret = 0.5*np.sin(xin-np.pi)*np.sin(yin)
    else:
        raise Exception('No such case, u2')

    return ret 

# RHS update
@njit
def update_inner(uin_inner,dusol):
    myu = uin_inner.copy()
    myu[buffer_x:Nx-buffer_x,buffer_y:Ny-buffer_y] =  \
        uin_inner[buffer_x:Nx-buffer_x,buffer_y:Ny-buffer_y]  \
          + dusol[buffer_x:Nx-buffer_x,buffer_y:Ny-buffer_y] 
    return myu

# Get averaged quantity over transverse diretion
@njit
def get_bar(uin):
    ubar = np.zeros(Nx)
    ubar_full = np.zeros_like(uin)
    for i in range(Nx):
        for j in range(Ny):
            ubar[i] += uin[i,j]
        ubar[i] = ubar[i]/Ny

    for i in range(Nx):
        for j in range(Ny):
            ubar_full[i,j] = ubar[i]

    return ubar, ubar_full


# Ops
@njit
def dudx(uin,i,j):
    ret = (uin[i+1][j] - uin[i-1][j])/(2.*dx)
    return ret

@njit
def dudy(uin,i,j):
    ret = (uin[i][j+1] - uin[i][j-1])/(2.*dy)
    return ret

@njit
def du2dx2(uin,i,j):
    ret = (uin[i+1][j] - 2*uin[i][j] + uin[i-1][j])/(dx*dx)
    return ret

@njit
def du2dy2(uin,i,j):
    ret = (uin[i][j+1] - 2*uin[i][j] + uin[i][j-1])/(dy*dy)
    return ret

@njit
def laplace2d(uin,i,j):
    ret = (uin[i+1][j] + uin[i-1][j] + \
           uin[i][j+1] + uin[i][j-1] - \
           4*uin[i][j]) / (dx*dy)
    return ret

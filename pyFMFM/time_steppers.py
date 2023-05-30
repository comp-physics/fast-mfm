""" 
Time steppers here
"""

import numpy as np
from ops import update_inner, rhs
from params import *
from numba import njit
    
@njit
def euler(uin):
    du = dt*rhs(uin)
    # uout = update_inner(uin,du)
    # return uout
    return du

@njit
def rk2(uin):
    k1 = dt*rhs(uin)
    k2 = dt*rhs(uin+k1)
    du = (k1+k2)/2.
    # uout = update_inner(uin,du)
    # return uout
    return du

@njit
def rk4(uin):
    k1 = dt*rhs(uin)
    k2 = dt*rhs(uin+k1/2.)
    k3 = dt*rhs(uin+k2/2.)
    k4 = dt*rhs(uin+k3)
    du = (k1+2*k2+2*k3+k4)/6.
    # uout = update_inner(uin,du)
    # return uout
    return du

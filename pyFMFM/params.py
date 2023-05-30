""" 
This is where the global variables go
"""

import numpy as np

# save_stuff = False
save_stuff = True

# make_plots = False
make_plots = True

make_animation = False
# make_animation = True

apply_mfm = False
# apply_mfm = True


## For running MFM with multiple columns of perturbation
perturb_cbars = False
# perturb_cbars = True


## For running MFM until sbar is steady
s_term_crit = False
# s_term_crit = 1e-9

# if s_term_crit and save_stuff:
#     raise Exception('Case error, save_stuff and s_term_crit')


# case = 'heat'
# case = 'advection'
# case = 'advection-diffusion'
# case = 'parallel'
# case = 'inhomogeneous-periodic'
case = 'inhomogeneous'

# tstepper = 'euler'
# tstepper = 'rk2'
tstepper = 'rk4'

bc_x = 'periodic'

if case == 'inhomogeneous':
    bc_y = 'noflux'
elif case == 'inhomogeneous-periodic':
    bc_y = 'periodic'
    # bc_y = 'noflux'
else:
    bc_y = 'none'

# ic = 'zeros'
ic = '1dgauss'
# ic = '2dgauss'
# ic = 'hotend'
# ic = 'florian'
# ic = 'loopy'

Nx = 100
Ny = 100

# Lx = 8.
# Ly = 4.

# Lx = 4.
Lx = 2.*np.pi
Ly = 2.*np.pi

x = np.linspace(0,Lx,Nx)
y = np.linspace(0,Ly,Ny)

xv, yv = np.meshgrid(x, y, indexing='ij')

dx = x[1]-x[0]
dy = y[1]-y[0]

adv_speed = 1.

cfl=0.1
dt = cfl*dx*dy/adv_speed

if s_term_crit:
    Nt = int(5*1e4)
    T = dt*Nt
else:
    T = 1.
    Nt = int(T/dt)
    # Nt = 100

Nt_stride = 100

myt = np.zeros(Nt)
for k in range(Nt):
    myt[k] = k*dt

t1t, x1t = np.meshgrid(myt, x, indexing='ij')

if bc_x == 'periodic' or bc_x == 'none' or bc_x == 'noflux':
    buffer_x = 1
else:
    raise Exception('That x-dir BC doesnt exist')

if bc_y == 'periodic' or bc_y == 'none' or bc_y == 'noflux':
    buffer_y = 1
else:
    raise Exception('That y-dir BC doesnt exist')

if ic == '1dgauss' or ic == '2dgauss':
    plotmax = 1
elif case == 'hotend':
    plotmax = 100
else:
    plotmax = 1

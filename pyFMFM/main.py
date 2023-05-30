#!/usr/bin/env python3

"""
Driver - call this as an executable via Python3
"""

import numpy as np
import time

from params import *
from time_steppers import *
from multicolor import *
from plot import *
from ic import *
from ops import *
from loop import *

import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation

import argparse
 
parser = argparse.ArgumentParser(description="Just an example",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-x", "--xwrite",  action="store_true", help="archive mode")
parser.add_argument("-r", "--rhsread", action="store_true", help="archive mode")
parser.add_argument('-l', "--loc", action='store', type=str, help='The text to parse.')
args = parser.parse_args()
config = vars(args)
locals().update(config)
print(config)

if xwrite:
    print('Writing x')
    np.savetxt(loc + 'x.csv', x[buffer_x:Nx-buffer_x], delimiter=',')
    raise Exception("Wrote csv of x")
elif rhsread:
    print('Reading IC')
    sources = np.loadtxt(loc + 'ic.csv', delimiter=',')
    ic = 'loopy'
    my_sbar = multi_input_timestep(ics = sources)
    np.savetxt(loc + 'sbar.csv', my_sbar, delimiter=',')
else:
    timestep_wrapper()
    # raise Exception("Nothing to do")


if make_animation and save_stuff:
    def animate(k):
        plotheatmap(u_save[k], k, save=False)

    file_name ='D/sol' 
    interv = 1000
    anim = animation.FuncAnimation(plt.figure(), animate, interval=interv, frames=range(0,Nt,interv), repeat=False)
    anim.save(file_name+'.gif')

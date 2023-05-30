"""
Main timestepping loop, saving of solution in all of time
"""

import numpy as np

from params import *
from time_steppers import *
from ic import *
from ops import *
from multicolor import *
from plot import *
import time


if tstepper == 'euler':
    stepper = euler
elif tstepper == 'rk2':
    stepper = rk2
elif tstepper == 'rk4':
    stepper = rk4
else:
    raise Exception('no such time stepper')

if save_stuff:
    u_save = np.zeros((Nt, Nx, Ny))

if apply_mfm and save_stuff:
    ubar_save    = np.zeros_like(u_save)
    uoldbar_save = np.zeros_like(u_save)
    unewbar_save = np.zeros_like(u_save)
    s_save       = np.zeros_like(u_save)
    dubar_save   = np.zeros_like(u_save)
    ds_save      = np.zeros_like(u_save)
    duxt_save    = np.zeros(Nt)
    dsxt_save    = np.zeros(Nt)

def multi_input_timestep(ics):
    Ninputs = ics.shape[1]
    sbar_out = np.zeros_like(ics) 

    xlen_expected = Nx-2*buffer_x
    if xlen_expected != ics.shape[0]:
        raise Exception('Wrong size IC from Florian')

    print('Nin =',Ninputs)
    for i in range(Ninputs):
        u = np.zeros((Nx, Ny))
        u = set_ic(u,coloring=ics[:,i])
        u, sbar = timestep(u) 
        sbar_out[:,i] = sbar.copy()

    return sbar_out

def timestep_wrapper():
    u = np.zeros((Nx, Ny))

    if apply_mfm and perturb_cbars:
        print('Running MFM with several perturbed ICs')
        Nradius = 3
        picked_column = 3
        ranges, centers = ranges_and_centers(Nradius, Nx)

        # Perturb serveral separated columns at once
        cbar = construct_cbar(centers, Nx)
        u = set_ic(u,coloring=cbar)
    else:
        print('Running with regular IC')
        u = set_ic(u)

    # plt.clf()
    # plt.xlabel("x")
    # plt.ylabel("y")
    # plt.pcolormesh(xv, yv, u, cmap=plt.cm.jet, vmin=0, vmax=plotmax)
    # plt.colorbar()
    # plt.savefig('florianos_full_ic.png')
    # plt.close()

    if apply_mfm:
        u, mysbar = timestep(u) 
        if perturb_cbars:
            sbars = deconstruct_sbar(ranges, mysbar)
    else:
        u = timestep(u)

    if False:
        # Pick just one center/column
        if apply_mfm:
            sbars_save = sbars.copy()
            u = np.zeros((Nx, Ny))
            ranges, centers = ranges_and_centers(Nradius, Nx)
            cbar = construct_cbar(centers, Nx, pick=picked_column)
            u = set_ic(u,coloring=cbar)

            # plt.clf()
            # plt.xlabel("x")
            # plt.ylabel("y")
            # plt.pcolormesh(xv, yv, u, cmap=plt.cm.jet, vmin=0, vmax=plotmax)
            # plt.colorbar()
            # plt.savefig('florianos_picked_ic.png')
            # plt.close()

            u, mysbar = timestep(u) 

            print('Picked sbars =', picked_column,  sbars_save[:,picked_column])
            print('Mysbar =',mysbar)
            # plot_centers = np.zeros((2,len(centers)))
            # plot_centers[0,:] = centers[:]
            # print('plot_centers',plot_centers)
            
            plt.clf()
            # plt.xlabel("t")
            # plt.ylabel("Linf(dubar[t])")
            plt.plot(mysbar)
            plt.plot(sbars_save[:,picked_column])
            plt.plot(centers,np.zeros(len(centers)),'o')
            plt.savefig('sbars.png')
            plt.close()

    # raise Exception('done')

def timestep(u):
    if apply_mfm:
        ubar = np.zeros(Nx)
        ubar_full = np.zeros((Nx, Ny))
        s = np.zeros_like(u)

    # plt.clf()
    # plt.xlabel("x")
    # plt.ylabel("y")
    # plt.pcolormesh(u, cmap=plt.cm.jet, vmin=0, vmax=plotmax)
    # plt.colorbar()
    # plt.savefig('ic.png')
    # plt.close()

    # Step through
    for k in range(Nt):
        if k%Nt_stride == 0:
            print('Step:',k,'of',Nt)

        if save_stuff:
            u_save[k,:,:] = u[:,:].copy() 

        # save original stuff
        if apply_mfm:
            if save_stuff:
                ubar_save[k,:,:] = ubar_full[:,:].copy() 
                s_save[k,:,:] = s[:,:].copy() 

            sold = s.copy()
            uold = u.copy()
            ubar, ubar_full = get_bar(u)
        
            if save_stuff:
                uoldbar_save[k,:,:] = ubar_full[:,:].copy() 

        # Gets the "effective" RHS
        du = stepper(u)
    
        # actually take the "time step", update u
        u = update_inner(u,du)

        if apply_mfm:
            # compute forcing
            s = u - uold

            # average forcing
            sbar, sbar_full = get_bar(s)
            # apply the forcing as control
            u = u - sbar_full

            # get updated ubar for comparison
            unewbar, unewbar_full = get_bar(u)

            # compute some per-time-step differences
            dubar = np.max(np.abs(unewbar[:] - ubar[:]))
            dsbar = np.max(np.abs(s - sold))

            if save_stuff:
                unewbar_save[k,:,:] = unewbar_full[:,:].copy() 
                dubar_save[k,:,:] = unewbar_save[k,:,:] - uoldbar_save[k,:,:]
                ds_save[k,:,:] = s[:,:] - sold[:,:]
                duxt_save[k] = dubar
                dsxt_save[k] = dsbar

            # print them
            if k%Nt_stride == 0:
                print('dsbar =', dsbar)
                # print('dubar =', dubar)

            if dsbar < s_term_crit:
                print('dsbar < scrit, terminate')
                break 

    if make_plots and save_stuff:
        plotheatmap(u_save[Nt-1],Nt-1,save=True)

    if make_plots and apply_mfm and save_stuff:
        plt.clf()
        plt.xlabel("t")
        plt.ylabel("Linf(dubar[t])")
        plt.plot(myt,duxt_save)
        plt.savefig('D/dubar-of-time.png')
        plt.close()

        plt.clf()
        plt.xlabel("t")
        plt.ylabel("Linf(dsbar[t])")
        plt.semilogy(myt,dsxt_save)
        plt.savefig('D/dsbar-of-time.png')
        plt.close()

        plt.clf()
        plt.xlabel("x")
        plt.ylabel("t")
        plt.title("sbar (log scale)")
        plt.pcolormesh(x1t, t1t, np.log10(np.abs(s_save[:,:,0])+1e-14), cmap=plt.cm.pink)
        plt.colorbar()
        plt.savefig('D/stime.png')
        plt.close()

        plt.clf()
        plt.xlabel("x")
        plt.ylabel("t")
        plt.pcolormesh(x1t, t1t, unewbar_save[:,:,0], cmap=plt.cm.pink)
        plt.colorbar()
        plt.title("ubar")
        plt.savefig('D/ubartime.png')
        plt.close()

        plt.clf()
        plt.xlabel("x")
        plt.ylabel("t")
        plt.title("dubar = ubar[t_i] - ubar[t_(i-1)]")
        plt.pcolormesh(x1t, t1t, dubar_save[:,:,0], cmap=plt.cm.pink)
        plt.colorbar()
        plt.savefig('D/dubartime.png')
        plt.close()


        plt.clf()
        plt.xlabel("x")
        plt.ylabel("t")
        plt.title("ds = s[t_i]-s[t_(i-1)] (log scale)")
        plt.pcolormesh(x1t, t1t, np.log10(np.abs(ds_save[:,:,0])+1e-14), cmap=plt.cm.pink)
        plt.colorbar()
        plt.savefig('D/dstime.png')
        plt.close()
    
    if apply_mfm:
        return u, sbar[buffer_x:Nx-buffer_x]
    else:
        return u

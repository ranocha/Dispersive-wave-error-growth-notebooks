#!/usr/bin/env python
# encoding: utf-8
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import os
from six.moves import range
from clawpack.petclaw.solution import Solution
from clawpack.petclaw.fileio.petsc import write
from matplotlib import pyplot as plt
import sys
sys.path.append('./../')
import psystem
import h5py

def qinit(state,ref_state):
    # data structures for initial condition
    x = state.grid.x.centers
    dx=x[1]-x[0]
    mx = len(x)
    q0 = np.zeros_like(state.q[0,:])
    q1 = np.zeros_like(state.q[1,:])
    
    # reference solution
    ref_x = ref_state.grid.x.centers
    ref_mx = len(ref_x)
    ref_q0 = ref_state.q[0,:]
    ref_q1 = ref_state.q[1,:]

    # do high-order reconstruction
    ref_q0m = np.zeros_like(ref_q0)
    ref_q1m = np.zeros_like(ref_q1)
    ref_q0m[2:-2] = (9*ref_q0[:-4]-116*ref_q0[1:-3]+2134*ref_q0[2:-2]-116*ref_q0[3:-1]+9*ref_q0[4:])/1920.0
    ref_q1m[2:-2] = (9*ref_q1[:-4]-116*ref_q1[1:-3]+2134*ref_q1[2:-2]-116*ref_q1[3:-1]+9*ref_q1[4:])/1920.0

    # construct interpolant
    from scipy.interpolate import interp1d
    from scipy.integrate import quad
    ref_q0Interp = interp1d(ref_x,ref_q0m)
    ref_q1Interp = interp1d(ref_x,ref_q1m)

    # zero order extraplation
    factor = int(mx/ref_mx)
    for i in range(factor):
        q0[i::factor] = ref_q0[:]
        q1[i::factor] = ref_q1[:]
    #

    # Higher-order interpolation
    # Get starting and ending points for integration
    #amax = q0.argmax()
    #tol = 0
    #for k in reversed(xrange(amax)):
    #    if q0[k] <= tol:
    #        kStart = k - 50 # with safety factor
    #        break
    #    #
    #
    #for k in xrange(amax,mx):
    #    if q0[k] <= tol:
    #        kEnd = k + 50 # with safety factor
    #        break
    #    #
    #

    # get averages per cell based on high-order interpolation
    #q0[kStart:kEnd] = np.array([1/dx*quad(ref_q0Interp,x[k]-dx/2.,x[k]+dx/2.)[0] for k in range(kStart,kEnd)])
    #q1[kStart:kEnd] = np.array([1/dx*quad(ref_q1Interp,x[k]-dx/2.,x[k]+dx/2.)[0] for k in range(kStart,kEnd)])

    # save refined solution in state
    state.q[0,:] = q0
    state.q[1,:] = q1
#

def setaux(x,
           rhoB=4,
           KB=4,
           rhoA=1,
           KA=1,
           alpha=0.5,
           coeff_type='smooth'):

    aux = np.empty([3,len(x)],order='F')
    xfrac = x-np.floor(x)

    if coeff_type == 'smooth':
        #Density:
        aux[0,:] = (rhoA+rhoB)/2 + (rhoA-rhoB)/2*np.sin(2*np.pi*x)
        #Bulk modulus:
        aux[1,:] = (KA+KB)/2 + (KA-KB)/2*np.sin(2*np.pi*x)
        aux[2,:] = 0. # not used
    else:
        #Density:
        aux[0,:] = rhoA*(xfrac<alpha)+rhoB*(xfrac>=alpha)
        #Bulk modulus:
        aux[1,:] = KA  *(xfrac<alpha)+KB  *(xfrac>=alpha)
        aux[2,:] = 0. # not used
    #
    return aux

def refine(outdir='./_output',
           file_prefix=None,
           frame=400,
           # Refinement
           refn=1,
           Lx=400):
    import clawpack.petclaw as pyclaw
    
    ref_sol=Solution(frame,file_format='petsc',read_aux=False,path='./_output_cut_steg',file_prefix='cut_steg')
    ref_xc = ref_sol.state.grid.x.centers
    
    # Domain and mesh
    xlower=0.0; xupper=Lx
    mx = len(ref_xc) * (2**refn)
    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    domain = pyclaw.Domain(x)

    # State
    num_eqn = ref_sol.state.num_eqn
    num_aux = ref_sol.state.num_aux
    state = pyclaw.State(domain,num_eqn,num_aux=num_aux)

    # Problem data
    t1 = ref_sol.state.problem_data['t1']
    tw1 = ref_sol.state.problem_data['tw1']
    amp_wall = ref_sol.state.problem_data['amp_wall']
    alpha = ref_sol.state.problem_data['alpha']
    KA = ref_sol.state.problem_data['KA']
    KB = ref_sol.state.problem_data['KB']
    rhoA = ref_sol.state.problem_data['rhoA']
    rhoB = ref_sol.state.problem_data['rhoB']
    coeff_type = ref_sol.state.problem_data['coeff_type']

    state.problem_data = {}
    state.problem_data['stress_relation'] = 'exponential'
    state.problem_data['t1']    = t1
    state.problem_data['tw1']   = tw1
    state.problem_data['amp_wall'] = amp_wall
    state.problem_data['alpha'] = alpha
    state.problem_data['KA'] = KA
    state.problem_data['KB'] = KB
    state.problem_data['rhoA'] = rhoA
    state.problem_data['rhoB'] = rhoB
    state.problem_data['coeff_type'] = coeff_type
    state.problem_data['trtime'] = 999999999.0
    state.problem_data['trdone'] = False

    #Initialize q and aux
    xc=state.grid.x.centers
    state.aux[:,:] = setaux(xc,rhoB,KB,rhoA,KA,alpha,coeff_type=coeff_type)

    # initial condition 
    qinit(state,ref_sol.state)

    # create solution
    sol = pyclaw.Solution(state,domain)

    # output 
    if not os.path.exists('./_output'): os.mkdir('_output')
    write(sol,0,'./_output/',file_prefix='init_refn'+str(refn))

    
    #######################################
    # ***** plot reference solution ***** #
    #######################################
    eps=ref_sol.state.q[0,:]
    sigma = psystem.stress(ref_sol.state)
    plt.figure(figsize=(15,5))
    xmax = ref_xc[sigma.argmax()]
    plt.plot(ref_xc,eps,'-r',lw=1)    
    plt.plot(ref_xc,sigma,'-k',lw=3)    
    plt.xlim([xmax-10,xmax+10])
    plt.ylim([-0.1,0.7])
    plt.savefig('refn0.png')

    #####################################
    # ***** plot refined solution ***** #
    #####################################
    eps=sol.state.q[0,:]
    vel=sol.state.q[1,:]
    sigma = psystem.stress(sol.state)
    plt.figure(figsize=(15,5))
    xmax = xc[sigma.argmax()]
    plt.plot(xc,eps,'-r',lw=1)    
    plt.plot(xc,sigma,'-k',lw=3)
    #plt.plot(ref_xc,ref_sol.state.q[1,:],'-b',lw=1)    
    plt.xlim([xmax-10,xmax+10])
    plt.ylim([-0.1,0.7])
    plt.savefig('refn'+str(refn)+'.png')

    #############################################
    # ***** save refined solution to hdf5 ***** #
    #############################################
    name='refn'+str(refn)
    # get left and right points
    amax = sigma.argmax()
    dx=xc[1]-xc[0]
    axL = np.abs(np.floor(xc[amax]-10)-(xc-dx/2.0)).argmin()
    axR = np.abs(np.floor(xc[amax]+10)-(xc-dx/2.0)).argmin()
    # shift location
    xShift = xc[axL:axR]-xc[axL]+dx/2.0
    hf = h5py.File(name+'_small_domain.h5', 'w')
    hf.create_dataset('x', data=xShift)
    hf.create_dataset('stress', data=sigma[axL:axR])
    hf.create_dataset('strain', data=eps[axL:axR])
    hf.create_dataset('vel', data=vel[axL:axR])
    hf.close()
    
    hf2 = h5py.File(name+'_full_domain.h5', 'w')
    hf2.create_dataset('x', data=xc)
    hf2.create_dataset('stress', data=sigma)
    hf2.create_dataset('strain', data=eps[:])
    hf2.create_dataset('vel', data=vel[:])
    hf2.close()
    #

    ##################################################################
    # ***** Write clawpack solution in the small domain (0,20) ***** #
    ##################################################################
    xlower=0.0; xupper=20.0
    mx = len(xShift)
    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    domain = pyclaw.Domain(x)
    
    # State
    state = pyclaw.State(domain,num_eqn,num_aux=num_aux)

    # Problem data
    state.problem_data = {}
    state.problem_data['stress_relation'] = 'exponential'
    state.problem_data['t1']    = t1
    state.problem_data['tw1']   = tw1
    state.problem_data['amp_wall'] = amp_wall
    state.problem_data['alpha'] = alpha
    state.problem_data['KA'] = KA
    state.problem_data['KB'] = KB
    state.problem_data['rhoA'] = rhoA
    state.problem_data['rhoB'] = rhoB
    state.problem_data['coeff_type'] = coeff_type
    state.problem_data['trtime'] = 999999999.0
    state.problem_data['trdone'] = False

    #Initialize q and aux
    xc=state.grid.x.centers
    state.aux[:,:] = setaux(xc,rhoB,KB,rhoA,KA,alpha,coeff_type=coeff_type)

    # initial condition
    state.q[0,:] = eps[axL:axR]
    state.q[1,:] = vel[axL:axR]

    # create solution
    sol = pyclaw.Solution(state,domain)

    # output 
    if not os.path.exists('./_output_small_domain'): os.mkdir('_output_small_domain')
    write(sol,0,'./_output_small_domain/',file_prefix='init_refn'+str(refn))
#

# reference solution: Nx=1024
refine(refn=0) # Nx=1024
refine(refn=1) # Nx=2048
refine(refn=2) # Nx=4096
refine(refn=3) # Nx=8192



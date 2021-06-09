#!/usr/bin/env python
# encoding: utf-8
r"""
Solitary wave formation in periodic nonlinear elastic media
===========================================================

Solve a one-dimensional nonlinear elasticity system:

.. math::
    \epsilon_t + u_x & = 0 \\
    (\rho(x) u)_t + \sigma(\epsilon,x)_x & = 0.

Here :math:`\epsilon` is the strain, :math:`\sigma` is the stress, 
u is the velocity, and :math:`\rho(x)` is the density.  
We take the stress-strain relation :math:`\sigma = e^{K(x)\epsilon}-1`;
:math:`K(x)` is the linearized bulk modulus.
Note that the density and bulk modulus may depend explicitly on x.

The problem solved here is based on [LeVYon03]_.  An initial hump
evolves into two trains of solitary waves.

"""
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
from six.moves import range


def qinit(state,wave_type,A=1.0,sigma2=10.0,xupper=600.):
    x = state.grid.x.centers
    
    if wave_type==0:
        # Gaussian
        sigma = A*np.exp(-(x-xupper/2.)**2/(2*sigma2))
        state.q[0,:] = np.log(sigma+1.)/state.aux[1,:]
        state.q[1,:] = 0.
    else:
        state.q[:,:] = 0.

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

    
def b4step(solver,state):
    #Reverse velocity at trtime
    #Note that trtime should be an output point
    if state.t>=state.problem_data['trtime']-1.e-10 and not state.problem_data['trdone']:
        #print 'Time reversing'
        state.q[1,:]=-state.q[1,:]
        state.q=state.q
        state.problem_data['trdone']=True
        if state.t>state.problem_data['trtime']:
            print('WARNING: trtime is '+str(state.problem_data['trtime'])+\
                ' but velocities reversed at time '+str(state.t))
    #Change to periodic BCs after initial pulse 
    if state.t>5*state.problem_data['tw1'] and solver.bc_lower[0]==0:
        solver.bc_lower[0] = 2
        solver.bc_upper[0] = 2
        solver.aux_bc_lower[0] = 2
        solver.aux_bc_upper[0] = 2

def zero_bc(state,dim,t,qbc,auxbc,num_ghost):
    """Set everything to zero"""
    if dim.on_upper_boundary:
        qbc[:,-num_ghost:]=0.

def moving_wall_bc(state,dim,t,qbc,auxbc,num_ghost):
    """Initial pulse generated at left boundary by prescribed motion"""
    if dim.on_lower_boundary:
        qbc[0,:num_ghost]=qbc[0,num_ghost] 
        t=state.t; t1=state.problem_data['t1']; tw1=state.problem_data['tw1']
        amp_wall=state.problem_data['amp_wall'];
        t0 = (t-t1)/tw1
        if abs(t0)<=1.: vwall = -amp_wall*(1.+np.cos(t0*np.pi))
        else: vwall=0.
        for ibc in range(num_ghost-1):
            qbc[1,num_ghost-ibc-1] = 2*vwall*state.aux[1,ibc] - qbc[1,num_ghost+ibc]

def p_system(use_petsc=True,
             kernel_language='Fortran',
             solver_type='classic',
             path='./_output',
             file_prefix=None,
             # Refinement
             refn=1, # 1=16
             # General parameters for simulatiuon #
             final_time=400.0,
             nDOut=400,
             restart_from_frame=None,
             # parameters for the coefficients
             rhoB=4.0,
             KB=4.0,
             rhoA=1.0,
             KA=1.0,
             alpha=0.5,
             # about the generation of the wave
             PERIODIC_BCs=False,
             name_file=None,
             wave_type = 0, # 0: via initial condition, 1: via left boundary
             amp_wall = 0.1,
             A=1.0,
             sigma2=10.0,
             coeff_type='smooth', 
             Lx=300,
             PERIODIC_DOMAIN=True):

    assert coeff_type in ['smooth','non-smooth']
    # ********** Import libraries ********** #
    from clawpack import riemann

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    # ********** Setup solver and solver parameters ********** #
    if kernel_language=='Python':
        rs = riemann.nonlinear_elasticity_1D_py.nonlinear_elasticity_1D
    elif kernel_language=='Fortran':
        rs = riemann.nonlinear_elasticity_fwave_1D

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D(rs)
        solver.char_decomp=0
        solver.limiters = pyclaw.limiters.tvd.minmod
    else:
        solver = pyclaw.ClawSolver1D(rs)
        solver.limiters = pyclaw.limiters.tvd.minmod
        #solver.limiters = pyclaw.limiters.tvd.superbee
        #solver.limiters = pyclaw.limiters.tvd.MC
    solver.kernel_language = kernel_language

    # ********** Boundary conditionos ********** #
    if wave_type == 0 or PERIODIC_BCs:
        solver.bc_lower[0] = pyclaw.BC.periodic
        solver.bc_upper[0] = pyclaw.BC.periodic
        
        #Use the same BCs for the aux array
        solver.aux_bc_lower[0] = pyclaw.BC.periodic
        solver.aux_bc_upper[0] = pyclaw.BC.periodic
    else:
        solver.bc_lower[0] = pyclaw.BC.custom
        solver.bc_upper[0] = pyclaw.BC.extrap
        
        #Use the same BCs for the aux array
        solver.aux_bc_lower[0] = pyclaw.BC.extrap
        solver.aux_bc_upper[0] = pyclaw.BC.extrap
    #

    claw = pyclaw.Controller()
    claw.tfinal = final_time
    claw.num_output_times = nDOut
    claw.solver = solver
    claw.outdir = path
    claw.keep_copy = True
    claw.output_style = 1
    claw.name_file = name_file
    
    solver.max_steps = 5000000
    solver.cfl_max = 0.11
    solver.cfl_desired = 0.1
    solver.fwave = True 
    solver.before_step = b4step 
    solver.user_bc_lower=moving_wall_bc
    solver.user_bc_upper=zero_bc
    
    # ********** Restart or create the simulation ********** #
    if restart_from_frame is not None:
        claw.solution = pyclaw.Solution(restart_from_frame, file_format='petsc',read_aux=False,file_prefix=file_prefix,path=path)
        grid = claw.solution.domain.grid
        xc = grid.x.centers
        claw.solution.state.aux[:,:] = setaux(xc,rhoB,KB,rhoA,KA,alpha,coeff_type=coeff_type)
        claw.num_output_times = claw.num_output_times - restart_from_frame
        claw.start_frame = restart_from_frame
    else:
        # Domain and mesh
        xlower=0.0; xupper=Lx
        Nx = 8 * (2**refn)
        mx=int((xupper-xlower)*Nx)
        x = pyclaw.Dimension(xlower,xupper,mx,name='x')
        domain = pyclaw.Domain(x)

        # State
        state = pyclaw.State(domain,solver.num_eqn,num_aux=3)        
    
        #Set global parameters
        state.problem_data = {}
        state.problem_data['stress_relation'] = 'exponential'
        state.problem_data['t1']    = 2.5
        state.problem_data['tw1']   = 2.5
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
        qinit(state,wave_type=wave_type,A=A,sigma2=sigma2,xupper=xupper)
                        
        claw.solution = pyclaw.Solution(state,domain)
    # run
    status = claw.run()
    #import pdb; pdb.set_trace()
    
def stress(state):
    """Compute stress from strain and momentum"""
    from clawpack.riemann.nonlinear_elasticity_1D_py import sigma
    # get data from state
    x = state.grid.x.centers
    rhoA = state.problem_data['rhoA']
    rhoB = state.problem_data['rhoB']
    KA = state.problem_data['KA']
    KB = state.problem_data['KB']
    alpha = state.problem_data['alpha']
    coeff_type = state.problem_data['coeff_type']

    # get coefficients
    aux = setaux(x,rhoB,KB,rhoA,KA,alpha,coeff_type=coeff_type)

    # get and return stress
    stress = sigma(state.q,aux,{'stress_relation':'exponential'})
    return stress


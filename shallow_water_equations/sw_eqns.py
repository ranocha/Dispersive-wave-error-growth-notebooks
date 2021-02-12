#!/usr/bin/env python
# encoding: utf-8

"""
Solve the 2D shallow water equations with bathymetry:
"""
import numpy as np

# ***************************************************************** #
# ********** GENERAL PARAMETERS USED FOR ALL SIMULATIONS ********** #
# ***************************************************************** #
# Boundary conditions and domain
Ly=1.0
PERIODIC_BCs = True

# Periodic properties of media
alphax=0.5; deltax=1000.0
alphay=0.5; deltay=1.0
b_A=0.0; b_B=0.5

# About the channel
bMax=0.5
yTilde=0.45 #yTilde\in[0,5] defines the inclination

def qinit(state,
          # about initial condition
          A,sig2,mwl,
          sw_diffracton=False):
    x0=0.
    y0=0.
    b = state.aux[0,:,:] # Bathymetry
    X,Y = state.grid.p_centers
    surface = mwl+A*np.exp(-X**2/(2*sig2))
    state.q[0,:,:] = surface - b
    state.q[1,:,:] = 0.
    state.q[2,:,:] = 0.

    if sw_diffracton:
        g = state.problem_data['grav']
        bm = 0.5*(b_A+b_B)
        bMwl=mwl-bm
        gamma=0
        x0=75.0
        absBeta=(b_A-b_B)**2/192.0/(mwl-b_A)/(mwl-b_B)
        sech = 1.0/np.cosh( np.sqrt(A/(8.0*absBeta/2.0*(1+gamma)*bMwl)) * (X-x0) )
        eta = A*sech**2 + mwl
        
        state.q[0,:,:] = eta - b
        state.q[1,:,:] = (eta-b)*( 2*np.sqrt(g*(eta-bm)) - 2*np.sqrt(g*bMwl) )
        state.q[2,:,:] = 0
    #
#
    
def bathymetry(x,y,coeff_type='smooth'):
    aux = np.empty((1,len(x),len(y)), order='F')

    if coeff_type == 'smooth':
        [yy,xx]=np.meshgrid(y,x)
        aux[0,:,:] = (b_A+b_B)/2.0 + (b_A-b_B)/2.0 * np.sin(2*np.pi*yy)
    else:
        # xfrac and yfrac are x and y relative to deltax and deltay resp.
        xfrac=x-np.floor(x/deltax)*deltax
        yfrac=y-np.floor(y/deltay)*deltay
        # create a meshgrid out of xfrac and yfrac
        [yyfrac,xxfrac]=np.meshgrid(yfrac,xfrac)
        # bathymetry
        aux[0,:,:] = (b_A*(xxfrac<=alphax*deltax)*(yyfrac<=alphay*deltay) +
                      b_A*(xxfrac >alphax*deltax)*(yyfrac >alphay*deltay) +
                      b_B*(xxfrac >alphax*deltax)*(yyfrac<=alphay*deltay) +
                      b_B*(xxfrac<=alphax*deltax)*(yyfrac >alphay*deltay))
    #
    return  aux
#

def switch_to_periodic_BCs(solver,state):
    from clawpack import pyclaw
    time_to_switch_BCs = state.problem_data['time_to_switch_BCs']

    #Change to periodic BCs after initial pulse
    if PERIODIC_BCs and state.t>time_to_switch_BCs and (solver.bc_lower[0]==pyclaw.BC.wall or solver.bc_lower[0]==pyclaw.BC.extrap):
        solver.bc_lower[0]=pyclaw.BC.periodic
        solver.bc_upper[0]=pyclaw.BC.periodic
        solver.aux_bc_lower[0]=pyclaw.BC.periodic
        solver.aux_bc_upper[0]=pyclaw.BC.periodic
#

def step_friction(solver, state, dt):
    "Friction source term:  -cf u / h.  This version is for Classic."
    cf = state.problem_data['cf']
    grav = state.problem_data['grav']
    q = state.q
    h = q[0,:,:]
    u = q[1,:,:]/h
    v = q[2,:,:]/h

    velNorm = np.sqrt(u**2+v**2)
    friction = grav*(cf**2)*velNorm/h**(4./3)

    q[1,:,:] = q[1,:,:] - dt*friction*h*u
    q[2,:,:] = q[2,:,:] - dt*friction*h*v
    #q[i,:,:] = q[i,:,:] - dt*cf*v/h
#

def sw_eqns(use_petsc=True,
            path='./_output',
            solver_type='classic',
            file_prefix=None,
            # Refinement
            refn=4, # goal: refn=4 -> Nx=128, Ny=128
            # General parameters for simulatiuon #
            final_time=400.0,
            nDOut=400,
            restart_from_frame=None,
            name_file=None,
            # about initial condition
            sw_diffracton=False,
            A=0.05,
            sig2=2,
            mwl=0.75,
            # about the bathymetry
            coeff_type='smooth', # smooth or non_smooth
            # about friction
            friction=False,
            friction_coeff=0.0,
            Lx=100.0,
            PERIODIC_DOMAIN = True,
            #switch to Periodic BCs
            time_to_switch_BCs = 25.0):
    #===========================================================================
    # Import libraries
    #===========================================================================
    from clawpack import riemann

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    #===========================================================================
    # Setup solver and solver parameters
    #===========================================================================
    if solver_type == 'classic':
        solver = pyclaw.ClawSolver2D(riemann.sw_aug_2D)
        #solver = pyclaw.ClawSolver2D(riemann.shallow_bathymetry_fwave_2D)
        solver.limiters = pyclaw.limiters.tvd.minmod
        solver.dimensional_split=True
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver2D(riemann.sw_aug_2D)
    #
    if sw_diffracton:
        solver.bc_lower[0] = pyclaw.BC.extrap
        #solver.bc_lower[0] = pyclaw.BC.wall
    else:
        solver.bc_lower[0] = pyclaw.BC.wall
    #
    solver.bc_upper[0] = pyclaw.BC.extrap
    solver.aux_bc_lower[0] = pyclaw.BC.extrap
    solver.aux_bc_upper[0] = pyclaw.BC.extrap

    if PERIODIC_DOMAIN==False:
        solver.bc_lower[1] = pyclaw.BC.wall
        solver.bc_upper[1] = pyclaw.BC.wall
        solver.aux_bc_lower[1] = pyclaw.BC.wall
        solver.aux_bc_upper[1] = pyclaw.BC.wall
    else:
        solver.bc_lower[1] = pyclaw.BC.periodic
        solver.bc_upper[1] = pyclaw.BC.periodic
        solver.aux_bc_lower[1] = pyclaw.BC.periodic
        solver.aux_bc_upper[1] = pyclaw.BC.periodic
    #

    solver.max_steps = 5000000
    solver.cfl_max = 0.11
    solver.cfl_desired = 0.1
    solver.fwave = True
    solver.before_step = switch_to_periodic_BCs

    #===========================================================================
    # Set up controller and controller parameters
    #===========================================================================
    claw = pyclaw.Controller()
    claw.tfinal = final_time
    claw.num_output_times = nDOut
    claw.solver = solver
    claw.outdir = path
    claw.keep_copy = True
    claw.name_file = name_file
 
    # Mesh resolution
    Nx = 8 * (2**refn)
    Ny = 8 * (2**refn)

    # FRICTION #
    if friction:
        solver.step_source = step_friction
        solver.source_split = 1
    #
    if restart_from_frame is not None:
        claw.solution = pyclaw.Solution(restart_from_frame, file_format='petsc',read_aux=False,file_prefix=file_prefix,path=path)
        grid = claw.solution.domain.grid
        claw.solution.state.aux = bathymetry(grid.x.centers,grid.y.centers,coeff_type=coeff_type)
        claw.num_output_times = claw.num_output_times - restart_from_frame
        claw.start_frame = restart_from_frame
        claw.solution.state.problem_data['time_to_switch_BCs'] = time_to_switch_BCs
    else:
        # Domain:
        xlower = 0.; xupper = Lx
        if PERIODIC_DOMAIN:
            ylower = -Ly/2.; yupper =  Ly/2.
        else:
            ylower = -Ly/4.; yupper =  Ly/4.
        #
        mx=int((xupper-xlower)*Nx)
        my=int((yupper-ylower)*Ny)
        
        x = pyclaw.Dimension(xlower,xupper,mx,name='x')
        y = pyclaw.Dimension(ylower,yupper,my,name='y')
        domain = pyclaw.Domain([x,y])

        num_aux = 1
        state = pyclaw.State(domain,solver.num_eqn,num_aux)
        state.aux = bathymetry(state.grid.x.centers,state.grid.y.centers,coeff_type=coeff_type)

        state.problem_data['grav'] = 9.8
        state.problem_data['cf'] = friction_coeff

        claw.solution = pyclaw.Solution(state,domain)
        qinit(state,A,sig2,mwl,sw_diffracton=sw_diffracton)
        state.problem_data['time_to_switch_BCs'] = time_to_switch_BCs

    # run
    status = claw.run()
#

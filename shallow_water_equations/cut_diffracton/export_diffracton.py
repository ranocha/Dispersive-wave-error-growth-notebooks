from clawpack.petclaw.solution import Solution
from clawpack.petclaw.fileio.petsc import write
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
from matplotlib import rc
import numpy as np
import h5py
import os

# ***************************************************************** #
# ********** GENERAL PARAMETERS USED FOR ALL SIMULATIONS ********** #
# ***************************************************************** #
# Periodic properties of media
alphax=0.5; deltax=1000.0
alphay=0.5; deltay=1.0
b_A=0.0; b_B=0.5

def qinit(state,ref_state):
    # data structures for initial condition
    x = state.grid.x.centers
    y = state.grid.y.centers
    mx = len(x)
    my = len(y)
    q0 = np.zeros_like(state.q[0,:,:])
    q1 = np.zeros_like(state.q[1,:,:])
    q2 = np.zeros_like(state.q[2,:,:])
    
    # reference solution
    ref_x = ref_state.grid.x.centers
    ref_y = ref_state.grid.y.centers
    ref_mx = len(ref_x)
    ref_my = len(ref_y)
    ref_q0 = ref_state.q[0,:,:]
    ref_q1 = ref_state.q[1,:,:]
    ref_q2 = ref_state.q[2,:,:]

    # refine by average
    x_factor = int(mx/ref_mx)
    y_factor = int(my/ref_my)
    for j in range(y_factor):
        for i in range(x_factor):
            q0[i::x_factor, j::y_factor] = ref_q0[:,:]
            q1[i::x_factor, j::y_factor] = ref_q1[:,:]
            q2[i::x_factor, j::y_factor] = ref_q2[:,:]
        #
    #
    # save refined solution in state
    state.q[0,:,:] = q0
    state.q[1,:,:] = q1
    state.q[2,:,:] = q2
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

def refine(
        frame = 340,
        # Refinement
        refn=0,
        # General parameters #
        coeff_type='smooth',
        mwl=0.75,
        Lx=100.0,
        plot_pcolor=True):
    print "*************************************************************"
    import clawpack.petclaw as pyclaw

    ref_sol=Solution(frame,file_format='petsc',read_aux=False,path='./_output_cut_wave',file_prefix='cut_wave')
    ref_xc = ref_sol.state.grid.x.centers
    ref_yc = ref_sol.state.grid.y.centers
    ref_mx = len(ref_xc)
    ref_my = len(ref_yc)
    ref_bath = bathymetry(ref_xc, ref_yc, coeff_type=coeff_type)[0]
    print "mx, my of reference solution: ", ref_mx, ref_my
    
    # Domain and mesh
    Ly = 1.0
    xlower = 0.; xupper = Lx
    ylower = -Ly/2.; yupper =  Ly/2.
    #
    mx = ref_mx * (2**refn)
    my = ref_my * (2**refn)
    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    y = pyclaw.Dimension(ylower,yupper,my,name='y')
    domain = pyclaw.Domain([x,y])
    print "mx, my of refined solution: ", mx, my
    
    # state
    num_eqn = 3
    num_aux = 1
    state = pyclaw.State(domain,num_eqn,num_aux)
    
    # problem data
    state.problem_data['grav'] = ref_sol.state.problem_data['grav']
    state.problem_data['cf'] = ref_sol.state.problem_data['cf']
    state.problem_data['time_to_switch_BCs'] = ref_sol.state.problem_data['time_to_switch_BCs']

    # bathymetry
    state.aux = bathymetry(state.grid.x.centers,state.grid.y.centers,coeff_type=coeff_type)

    print ("refine solution")
    # initial condition
    qinit(state,ref_sol.state)

    # create solution
    sol = pyclaw.Solution(state,domain)

    #############################################################
    # ***** refined solution vectors in the small doomain ***** #
    #############################################################
    xc = state.grid.x.centers
    yc = state.grid.y.centers
    dx=xc[1]-xc[0]
    # get the location of the max
    eta = state.q[0,:,:] + state.aux[0]
    amax = eta[:,my/2].argmax()
    # get left and right points
    axL = np.abs(np.floor(xc[amax]-10)-(xc-dx/2.0)).argmin()
    axR = np.abs(np.floor(xc[amax]+10)-(xc-dx/2.0)).argmin()
    # get x shifted
    xShifted = xc[axL:axR] - xc[axL]+dx/2.0
    # refined solution vector in the small domain
    h = state.q[0,axL:axR,:]
    hu = state.q[1,axL:axR,:]
    hv = state.q[2,axL:axR,:]

    #######################################
    # ***** plot reference solution ***** #
    #######################################
    ref_eta = ref_sol.state.q[0,:,:] + ref_bath
    print ("plot reference solution in the original domain")
    if plot_pcolor:
        (ref_yy,ref_xx) = np.meshgrid(ref_yc,ref_xc)
        pl.figure(figsize=(15,5))        
        pl.pcolormesh(ref_xx,ref_yy,ref_eta,cmap=pl.get_cmap('Blues'))
        pl.title("t="+str(ref_sol.state.t),fontsize=25)
        pl.xticks(size=25); pl.yticks(size=25)
        pl.xlim([75,90])
        cb = pl.colorbar()
        imaxes = pl.gca(); pl.axes(cb.ax)
        pl.yticks(fontsize=25); pl.axes(imaxes)
        pl.savefig("ref_soln.png")
        pl.close()
    
    # plot slices
    pl.figure(figsize=(15,5))
    mean_ref_eta=np.sum(ref_eta,1)/ref_my
    pl.plot(ref_xc,mean_ref_eta,'-k',lw=3)
    pl.plot(ref_xc,ref_eta[:,3*ref_my/4],'-r',lw=3)
    pl.plot(ref_xc,ref_eta[:,ref_my/4],'--b',lw=3)
    pl.title("t= "+str(ref_sol.state.t),fontsize=25)
    pl.xticks(size=25); pl.yticks(size=25)
    pl.tight_layout()
    pl.xlim([75,90])
    pl.gca().ticklabel_format(useOffset=False)
    pl.savefig("ref_soln_slices.png")
    pl.close()

    #####################################
    # ***** plot refined solution ***** #
    #####################################
    print ("plot refined solution in the smaller domain")
    if plot_pcolor:
        (yy,xx) = np.meshgrid(yc,xShifted)
        pl.figure(figsize=(15,5))
        pl.pcolormesh(xx,yy,eta[axL:axR],cmap=pl.get_cmap('Blues'))
        pl.title("t="+str(state.t),fontsize=25)
        pl.xticks(size=25); pl.yticks(size=25)
        cb = pl.colorbar()
        imaxes = pl.gca(); pl.axes(cb.ax)
        pl.yticks(fontsize=25); pl.axes(imaxes)
        pl.savefig("refn"+str(refn)+".png")
        pl.close()
    
    # plot slices
    pl.figure(figsize=(15,5))
    mean_eta=np.sum(eta[axL:axR],1)/my
    pl.plot(xShifted,mean_eta,'-k',lw=3)
    pl.plot(xShifted,eta[axL:axR,3*my/4],'-r',lw=3)
    pl.plot(xShifted,eta[axL:axR,my/4],'--b',lw=3)
    pl.title("t= "+str(sol.state.t),fontsize=25)
    pl.xticks(size=25); pl.yticks(size=25)
    pl.tight_layout()
    pl.gca().ticklabel_format(useOffset=False)
    pl.savefig("refn"+str(refn)+"_slices.png")
    pl.close()    
    
    #############################################
    # ***** save refined solution to hdf5 ***** #
    #############################################
    print ("write the refined solution in h5 files")
    name='refn'+str(refn)
    # shift location
    #xShift = xc[axL:axR]-xc[axL]+dx/2.0
    hf = h5py.File(name+'_small_domain.h5', 'w')
    hf.create_dataset('x', data=xShifted)
    hf.create_dataset('y', data=yc)
    hf.create_dataset('eta', data=eta[axL:axR,:])
    hf.create_dataset('h',  data=h)
    hf.create_dataset('hu', data=hu)
    hf.create_dataset('hv', data=hv)
    hf.close()
    
    ##################################################################
    # ***** Write clawpack solution in the small domain (0,20) ***** #
    ##################################################################
    print ("write refined solution in the small domain")
    
    # Domain and mesh
    Ly = 1.0
    xlower = 0.; xupper = 20.0
    ylower = -Ly/2.; yupper =  Ly/2.
    #
    mx = len(xShifted)
    my = len(yc) # same as for full domain
    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    y = pyclaw.Dimension(ylower,yupper,my,name='y')
    domain = pyclaw.Domain([x,y])
    
    # state
    state = pyclaw.State(domain,num_eqn,num_aux)
    
    # problem data
    state.problem_data['grav'] = ref_sol.state.problem_data['grav']
    state.problem_data['cf'] = ref_sol.state.problem_data['cf']
    state.problem_data['time_to_switch_BCs'] = ref_sol.state.problem_data['time_to_switch_BCs']

    # bathymetry
    state.aux = bathymetry(state.grid.x.centers,state.grid.y.centers,coeff_type=coeff_type)

    # initial condition
    state.q[0,:,:] = h
    state.q[1,:,:] = hu
    state.q[2,:,:] = hv

    # create solution
    sol = pyclaw.Solution(state,domain)

    # output
    print len(xShifted), len(yc)
    if not os.path.exists('./_output_small_domain'): os.mkdir('_output_small_domain')
    write(sol,0,'./_output_small_domain/',file_prefix='init_refn'+str(refn))
#

refine(refn=0) #Nx=Ny=512



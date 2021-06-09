import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
import numpy as np
import sys
sys.path.append('./../')
import psystem
from clawpack.petclaw.solution import Solution
from clawpack.petclaw.fileio.petsc import write

# tolerance for cutting the wave
tol=1.0E-12

def get_mask(sigma,x):
    mask = x*0

    mx=len(x)
    index_max=sigma.argmax()
    
    for k in reversed(xrange(index_max)):
        if (sigma[k] >= tol): 
            mask[k] = 1
        else:
            break

    for k in xrange(index_max,mx):
        if (sigma[k] >= tol): 
            mask[k] = 1
        else:
            break
        
    return mask
#

def cut(sol,mask):
    num_eqn = sol.state.num_eqn
    for i in xrange(num_eqn):
        sol.state.q[i,:]=sol.state.q[i,:]*mask
#

def plot(frame,
         pre_name='sigma_',
         file_prefix='claw',
         path='./_output/',
         plot_strain=True,
         xlimits=None,
         ylimits=None,
         follow=False,
         follow_window=10,
         mask=None):
    import sys
    sys.path.append('.')

    sol=Solution(frame,file_format='petsc',read_aux=False,path=path,file_prefix=file_prefix)
    x=sol.state.grid.x.centers
    eps=sol.state.q[0,:]

    # get stress
    sigma = psystem.stress(sol.state)

    # create figure and plot solution
    pl.figure(figsize=(15,5))
    if plot_strain:
        pl.plot(x,eps,'-r',lw=2)
    #
    pl.plot(x,sigma,'-k',lw=3)

    # To define the number for the name of the file
    if frame < 10:
        str_frame = "000"+str(frame)
    elif frame < 100:
        str_frame = "00"+str(frame)
    elif frame < 1000:
        str_frame = "0"+str(frame)
    else:
        str_frame = str(frame)
    #
    if follow:
        amax = sigma.argmax()
        xmax = x[amax]
        xlimits = [0,0]
        xlimits[0] = xmax-follow_window
        xlimits[1] = xmax+follow_window        
    #
    pl.title("t= "+str(sol.state.t),fontsize=25)
    pl.ylabel('$\sigma$',fontsize=25)
    pl.xticks(size=25); pl.yticks(size=25)
    if xlimits is not None:
        pl.xlim([xlimits[0],xlimits[1]])
    #
    if ylimits is not None:
        pl.ylim([ylimits[0],ylimits[1]])
    #
    if mask is not None:
        pl.plot(x,mask,'--b')
    #
    pl.savefig('./_plots/'+pre_name+'_'+str_frame+'.png')
    pl.close()
#

def cut_wave(path_solution='./_output/',
             file_prefix_solution='claw'):
    
    print 'reading the original solution'
    sol=Solution(frame,file_format='petsc',read_aux=False,path=path_solution,file_prefix=file_prefix_solution)
    x = sol.state.grid.x.centers
    eps=sol.state.q[0,:]
    sigma = psystem.stress(sol.state)
    
    print '... getting the mask'
    mask = get_mask(sigma,x)
    
    print '... cutting and writting the solitary wave'
    cut(sol,mask)
    write(sol,frame,'./_output_cut_steg/',file_prefix='cut_steg')
    
    print '... plotting the cut steg'
    plot(frame=frame,path='_output_cut_steg',file_prefix='cut_steg',pre_name='cut_steg')
    plot(frame=frame,path='_output_cut_steg',file_prefix='cut_steg',pre_name='cut_steg_follow',follow=True,follow_window=10,mask=mask)
    
#

if __name__== "__main__":
    import os
    if not os.path.exists('./_output_cut_steg/'): os.mkdir('./_output_cut_steg/')
    if not os.path.exists('./_plots/'): os.mkdir('./_plots/')

    import os
    frame=400
    cut_wave(path_solution='./_output/',
             file_prefix_solution='claw')

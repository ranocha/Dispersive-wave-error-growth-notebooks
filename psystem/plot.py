from clawpack.petclaw.solution import Solution
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
from matplotlib import rc
import numpy as np
import os

def plot_q(frame,
           file_prefix='claw',
           path='./_output/',
           xShift=0.0,
           xlimits=None,
           ylimits=None,
           name=None,
           plot_strain=True,
           plot_path='',
           ylabel='$\sigma$',
           plot_title=True,
           follow=False,
           follow_window=50,
           X=None, Eps=None, Vel=None, Sigma=None):
    import sys
    sys.path.append('.')
    import psystem

    sol=Solution(frame,file_format='petsc',read_aux=False,path=path,file_prefix=file_prefix)
    x=sol.state.grid.x.centers
    eps=sol.state.q[0,:]

    # get stress
    sigma = psystem.stress(sol.state)

    # number the frame to name the file
    if frame < 10:
        str_frame = "000"+str(frame)
    elif frame < 100:
        str_frame = "00"+str(frame)
    elif frame < 1000:
        str_frame = "0"+str(frame)
    else:
        str_frame = str(frame)

    # create the figure and plot the solution
    pl.figure(figsize=(15,5))
    if plot_strain:
        pl.plot(x+xShift,eps,'-r',lw=1)    
    pl.plot(x+xShift,sigma,'-k',lw=3)
    
    # format the plot
    if plot_title:
        pl.title("t= "+str(sol.state.t),fontsize=25)
    #
    if ylabel is not None:
        pl.ylabel(ylabel,fontsize=30)
    #
    pl.xticks(size=25); pl.yticks(size=25)
    if follow is True:
        amax = sigma.argmax()
        xmax = x[amax]
        xlimits = [0,0]
        xlimits[0] = xmax-follow_window
        xlimits[1] = xmax+follow_window
    if xlimits is not None:
        xlim=[xlimits[0], xlimits[1]]
    else:
        xlim=[np.min(x),np.max(x)]
    if ylimits is not None:
        ylim=[ylimits[0], ylimits[1]]
    else:
        ylim=[np.min(sigma),np.max(sigma)]
    pl.tight_layout()
    pl.axis([xlim[0]+xShift,xlim[1]+xShift,ylim[0],ylim[1]])
    pl.gca().ticklabel_format(useOffset=False)
    if name is None:
        if follow:            
            pl.savefig('./_plots_follow'+plot_path+'/sigma_'+str_frame+'.png',bbox_inches="tight")
        else:
            pl.savefig('./_plots'+plot_path+'/sigma_'+str_frame+'.png',bbox_inches="tight")
    else:
        pl.savefig(name+'.png',bbox_inches="tight")
    #
    pl.close()

    # save the data
    if X is not None:
        X.append(x)
        Eps.append(eps)
        Vel.append(sol.state.q[1,:])
        Sigma.append(sigma)
    #
    return x,eps,sigma
#

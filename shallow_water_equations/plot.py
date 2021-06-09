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
           plot_pcolor=True,
           plot_slices=True,
           plot_averaged_slices=False,
           xShift=0.0,
           slices_xlimits=None,
           slices_ylimits=None,
           pcolor_ylimits=None,
           name=None,
           plot_path='',
           ylabel='$\eta$',
           plot_title=True,
           legend=None,
           xKdV=None, uKdV=None, shiftKdV=0,
           X=None, eta1=None, eta2=None):
    import sys
    sys.path.append('.')
    import sw_eqns

    sol=Solution(frame,file_format='petsc',read_aux=False,path=path,file_prefix=file_prefix)
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers
    mx=len(x); my=len(y)

    h=sol.state.q[0,:,:]
    b=sw_eqns.bathymetry(x,y)[0,:,:]
    eta=h+b

    yy,xx = np.meshgrid(y,x)

    if frame < 10:
        str_frame = "000"+str(frame)
    elif frame < 100:
        str_frame = "00"+str(frame)
    elif frame < 1000:
        str_frame = "0"+str(frame)
    else:
        str_frame = str(frame)

    if plot_pcolor:
        pl.figure(figsize=(15,5))
        if slices_xlimits is not None:
            xlim=[slices_xlimits[0], slices_xlimits[1]]
        else:
            xlim=[np.min(x),np.max(x)]
        #
        if pcolor_ylimits is not None:
            ylim=[pcolor_ylimits[0], pcolor_ylimits[1]]
        else:
            ylim=[np.min(y), np.max(y)]
        #
        xdummy=np.linspace(xlim[0],xlim[1])
        #
        pl.pcolormesh(xx+xShift,yy,eta,cmap=pl.get_cmap('Blues'))
        if plot_title:
            pl.title("t= "+str(sol.state.t),fontsize=25)
        #
        pl.xticks(size=25); pl.yticks(size=25)
        cb = pl.colorbar()
        if slices_ylimits is not None:
            pl.clim(slices_ylimits[0],slices_ylimits[1])
        imaxes = pl.gca(); pl.axes(cb.ax)
        pl.yticks(fontsize=25); pl.axes(imaxes)
        pl.axis([xlim[0]+xShift,xlim[1]+xShift,ylim[0],ylim[1]])
        if name is None:
            pl.savefig('./_plots'+plot_path+'/eta_'+str_frame+'.png',bbox_inches="tight")
        else:
            pl.savefig(name+'.png',bbox_inches="tight")
        #
        pl.close()
    if plot_slices:
        pl.figure(figsize=(15,5))
        if plot_averaged_slices==True:
            mean_eta=np.sum(eta,1)/my
            #pl.plot(x+xShift,mean_eta,'-k',lw=3)
            x=np.concatenate([x,x[:25*128]+100])
            u=np.concatenate([mean_eta,mean_eta[:25*128]])
            pl.plot(x+xShift,u,'-k',lw=3)
            
        else:
            pl.plot(x+xShift,eta[:,3*my/4],'-r',lw=3)
            pl.plot(x+xShift,eta[:,my/4],'--b',lw=3)
            #
        #
        if plot_title:
            pl.title("t= "+str(sol.state.t),fontsize=25)
        #
        if ylabel is not None:
            pl.ylabel(ylabel,fontsize=30)
        #
        pl.xticks(size=25); pl.yticks(size=25)
        if slices_xlimits is not None:
            xlim=[slices_xlimits[0], slices_xlimits[1]]
        else:
            xlim=[np.min(x),np.max(x)]
        if slices_ylimits is not None:
            ylim=[slices_ylimits[0], slices_ylimits[1]]
        else:
            ylim=[np.min(eta),np.max(eta)]
        pl.tight_layout()
        pl.axis([xlim[0]+xShift,xlim[1]+xShift,ylim[0],ylim[1]])
        if uKdV is not None and xKdV is not None:
            pl.plot(xKdV+shiftKdV,uKdV[frame,:],'-r',lw=3)
            x=np.concatenate([xKdV,xKdV[:60]+200])
            u=np.concatenate([uKdV[125,:],uKdV[125,:60]])
            xdummy=np.linspace(300,310,100)
            #pl.plot(xdummy,xdummy*0+(np.sum(eta,1)/my)[-1],'-k',lw=3)
            #pl.plot(x+shiftKdV,u,'-r',lw=3)
        #
        if legend is not None:
            pl.legend(legend)
        #
        pl.gca().ticklabel_format(useOffset=False)
        if name is None:
            pl.savefig('./_plots'+plot_path+'/eta_'+str_frame+'_slices.png',bbox_inches="tight")
        else:
            pl.savefig(name+'_slices.png',bbox_inches="tight")
        #
        pl.close()

        if X is not None:
            X.append(x)
            eta1.append(eta[:,my/4])
            eta2.append(eta[:,3*my/4])
        #
    return eta
#

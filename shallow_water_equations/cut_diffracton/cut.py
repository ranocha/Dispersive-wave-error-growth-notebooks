from clawpack.petclaw.solution import Solution
from clawpack.petclaw.fileio.petsc import write
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl

import numpy as np

def bathymetry(x,y,b_A=0.,b_B=0.5,coeff_type='smooth'):
    aux = np.empty((1,len(x),len(y)), order='F')

    if coeff_type == 'smooth':
        [yy,xx]=np.meshgrid(y,x)
        aux[0,:,:] = (b_A+b_B)/2.0 + (b_A-b_B)/2.0 * np.sin(2*np.pi*yy)
    else:
        deltax=1000
        deltay=1.0
        alphax=0.5
        alphay=0.5
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

def get_mask(eta,x,y,mwl, plot_mask=True):
    yy,xx = np.meshgrid(y,x)
    mask = xx*0
    mask_inv = xx*0+1

    mx=len(x)
    my=len(y)

    tol = 1.0E-12
    index_max=eta[:,3*my/4].argmax()

    aux = eta - mwl
    print ('    ... reversed from the peak')
    for k in reversed(xrange(index_max)):
        #print (x[k],aux[k,3*my/4])
        if (aux[k,3*my/4] >= tol): 
            mask[k,:] = 1
            mask_inv[k,:] = 0
        else:
            break

    print ('    ... forward from the peak')        
    for k in xrange(index_max,mx):
        #print (x[k],aux[k,3*my/4])
        if (aux[k,3*my/4] >= tol): 
            mask[k,:] = 1
            mask_inv[k,:] = 0
        else:
            break

    if plot_mask:
        print ('    ... plotting mask')
        pl.pcolormesh(xx,yy,mask)
        pl.colorbar()
        pl.savefig('mask.png')

    return mask,mask_inv

def cut(sol,mwl,bathymetry,mask,mask_inv):
    num_eqn = sol.state.num_eqn
    for i in xrange(num_eqn):
        sol.state.q[i,:,:]=sol.state.q[i,:,:]*mask + (i==0)*(mwl-bathymetry)*mask_inv
#

def plot(frame,pre_name='eta_',file_prefix='claw',path='./_output/',slices_xlimits=None):
    import sys
    sys.path.append('.')

    sol=Solution(frame,file_format='petsc',read_aux=False,path=path,file_prefix=file_prefix)
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers
    mx=len(x); my=len(y)

    h=sol.state.q[0,:,:]
    b=bathymetry(x,y)[0,:,:]
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

    pl.figure(figsize=(15,5))
    pl.plot(x,eta[:,3*my/4],'-r',lw=3)
    pl.plot(x,eta[:,my/4],'--b',lw=3)
    pl.title("t= "+str(sol.state.t),fontsize=25)
    pl.xlabel('x',fontsize=25)
    pl.ylabel('$\eta$',fontsize=25)
    pl.xticks(size=25); pl.yticks(size=25)
    if slices_xlimits is not None:
        pl.axis([slices_xlimits[0],slices_xlimits[1],np.min(eta),np.max(eta)])
    pl.savefig('./_plots/'+pre_name+'_'+str_frame+'_slices.png')
    pl.close()

def cut_wave(mwl=0.75,
             path_solution='./_output/',
             file_prefix_solution='claw'):
    
    print '***** Cutting wave **********'
    print 'reading the original solution'
    
    sol=Solution(frame,file_format='petsc',read_aux=False,path=path_solution,file_prefix=file_prefix_solution)
    x = sol.state.grid.x.centers
    y = sol.state.grid.y.centers
    b = bathymetry(x,y)[0,:,:]
    eta = sol.state.q[0,:,:] + b

    my=len(y)
    pl.plot(x,eta[:,3*my/4])
    pl.savefig('aux.png')
    pl.close()
    
    print '... getting the mask'
    mask,mask_inv = get_mask(eta,x,y,mwl)
    
    print '... cutting and writting the solitary wave'
    cut(sol,mwl,b,mask,mask_inv)
    write(sol,frame,'./_output_cut_wave/',file_prefix='cut_wave')
    
    print '... plotting the isolated wave'
    plot(frame=frame,path='_output_cut_wave',file_prefix='cut_wave',pre_name='cut_wave')
    plot(frame=frame,path='_output_cut_wave',file_prefix='cut_wave',pre_name='cut_wave_follow',slices_xlimits=[70,90])
#

if __name__== "__main__":
    import os
    if not os.path.exists('./_output_cut_wave/'): os.mkdir('./_output_cut_wave/')
    if not os.path.exists('./_plots/'): os.mkdir('./_plots/')

    frame=340
    
    # cut wave 1
    cut_wave(mwl=0.75,
             path_solution='./_output/',
             file_prefix_solution='claw')


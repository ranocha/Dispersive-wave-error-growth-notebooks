#!/usr/bin/env python
# encoding: utf-8
import sys
sys.path.append('./../')

if __name__=="__main__":
    from plot import *

    for follow in [False,True]:
    
        if follow:
            if not os.path.exists('./_plots_follow'): os.mkdir('_plots_follow')
        else:
            if not os.path.exists('./_plots'): os.mkdir('_plots')
        #
        from_frame = 0
        to_frame   = 400
        frames=xrange(from_frame,to_frame+1)

        # data structures to save data
        X=[]
        Eps=[]
        Vel=[]
        Sigma=[]
        
        if not os.path.exists('./_plots'): os.mkdir('./_plots')
        print('**********************')
        print('**********************')
        print('Plotting solution ...')
        for i in frames:
            plot_q(frame=i,
                   path='./_output/',
                   xlimits=None,
                   ylimits=None, #[-0.02,0.625],
                   plot_strain=True,
                   follow=follow,
                   follow_window=50, #10, 50
                   X=X,
                   Eps=Eps,
                   Vel=Vel,
                   Sigma=Sigma)
            print ('frame '+str(i)+' plotted')
        #
        
    # Save solution to csv files
    np.savetxt("X_FV.csv",np.array(X),delimiter=',')
    np.savetxt("Eps_FV.csv",np.array(Eps),delimiter=',')
    np.savetxt("Vel_FV.csv",np.array(Vel),delimiter=',')
    np.savetxt("Sigma_FV.csv",np.array(Sigma),delimiter=',')


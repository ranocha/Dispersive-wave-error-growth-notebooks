#!/usr/bin/env python
# encoding: utf-8
import sys
sys.path.append('./../')

if __name__=="__main__":
    from sw_eqns import *
    sig2=2 #2, 10
    A=0.05
    
    sw_eqns(# General parameters for simulatiuon #
        final_time=340.0,
        nDOut=340,
        restart_from_frame=None,
        # About the refinement level
        refn=6, #6=512
        # About initial condition
        A=A,
        sig2=sig2,
        mwl=0.75)

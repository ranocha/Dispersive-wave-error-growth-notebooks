#!/usr/bin/env python
# encoding: utf-8
import sys
sys.path.append('./../')

if __name__=="__main__":
    from sw_eqns import *

    refn = 0
    path = './_output_refn' + str(refn)
        
    sw_eqns(# General parameters for simulatiuon #
        refn=refn,
        path=path,
        final_time=25.0,
        nDOut=25,
        restart_from_frame=0)

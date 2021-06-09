#!/usr/bin/env python
# encoding: utf-8
import sys
import os
sys.path.append('./../')
from petsc4py import PETSc

if __name__=="__main__":
    from psystem import *

    refn = 7 # goal 7 (1024 DoFs)
    outdir = './_output'
    rank = PETSc.COMM_WORLD.rank
    if rank==0:
        if not os.path.exists(outdir): os.mkdir(outdir)
    #
    p_system(refn=refn, #refn=1 -> Nx=16 DoFs per unit
             path=outdir,
             #solver_type='sharpclaw',
             final_time=400,
             nDOut=400,
             restart_from_frame=None,
             wave_type=1,
             A=1.0,
             sigma2=10.0,
             coeff_type='smooth',
             Lx=400)

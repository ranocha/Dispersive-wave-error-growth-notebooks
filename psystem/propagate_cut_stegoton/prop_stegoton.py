#!/usr/bin/env python
# encoding: utf-8
import sys
import os
sys.path.append('./../')
from petsc4py import PETSc

if __name__=="__main__":
    from psystem import *

    refn = 3
    path = './_output_refn' + str(refn)
    p_system(refn=refn, #refn=3 -> Nx=8192 DoFs per unit
             path=path,
             PERIODIC_BCs=True,
             final_time=25,
             nDOut=25,
             restart_from_frame=0,
             coeff_type='smooth')

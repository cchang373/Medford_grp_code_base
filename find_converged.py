# -*- coding: utf-8 -*-
"""
Created on Wed May 10 10:43:00 2017

@author: benjamin
"""
import os
from subprocess import call
import sys
sys.path.insert(0,"/nv/hp13/bcomer3/shared/espresso_rutile/tools")
from If_queue_allows import If_queue_allows

for root, dirs, files in os.walk(".", topdown=False):
    for name in files:
        TF = name == "converged_slab.traj"
        if TF:
            print(os.path.join(name,root))
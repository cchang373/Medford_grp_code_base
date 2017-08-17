# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 08:50:27 2017

@author: benjamin
"""
import os
from subprocess import call
#from subprocess import call
cur_dir=os.path.relpath(".","..")
import os
for root, dirs, files in os.walk(".", topdown=False):
    for name in files:
        if name == "converged_slab.traj":
            print(os.path.join(root, name))

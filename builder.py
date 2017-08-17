# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 09:01:59 2017

@author: benjamin

Cleans all restart.gpw files out of a given directory
"""
import os
from subprocess import call
#from subprocess import call
cur_dir=os.path.realpath(".")
import os
for root, dirs, files in os.walk(".", topdown=False):
    for name in files:
        if name == "build.py":
            print(os.path.join(root, name))
            os.chdir(os.path.join(root))
            call(["python","build.py"])
            os.chdir(cur_dir)

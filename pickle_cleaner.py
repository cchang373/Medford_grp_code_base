# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 09:01:59 2017

@author: benjamin

cleans empty vibrational pickles and submits vibration jobs that are either not
done or have never been started. Not rigorous
"""
import os
from pickle import load
from subprocess import call
import sys
sys.path.insert(0,"/nv/hp13/bcomer3/shared/espresso_rutile/tools")
from If_queue_allows import If_queue_allows
def pickle_cleaner():
    for root, dirs, files in os.walk(".", topdown=False):
        for name in files:
                if name.rsplit('.',1)[-1]=='pckl':
                          if os.stat(os.path.join(root,name)).st_size==0:
                              print os.path.join(root,name)
                              os.remove(os.path.join(root,name))


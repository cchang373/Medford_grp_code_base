# -*- coding: utf-8 -*-
"""
Created on Mon May  1 14:04:24 2017

@author: benjamin
"""

import os
from subprocess import call
import sys
sys.path.insert(0,"/nv/hp13/bcomer3/shared/espresso_rutile/tools")
import shutil


cur_dir=os.path.realpath(".")
import os
for root, dirs, files in os.walk(".", topdown=False):
    for name in files:
        TF = name == "qn_opt.py"
        if TF:
            shutil.copy('/home/benjamin/espresso/espresso_rutile/tools/qn_opt.py',os.path.join(root,name))
            
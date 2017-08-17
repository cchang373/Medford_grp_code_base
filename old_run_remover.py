# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 10:03:00 2017

@author: benjamin
"""

import os
from subprocess import call
cur_dir=os.path.relpath(".","..")
import os

#call(["source",  "/gpfs/pace1/project/chbe-medford/medford-share/envs/espresso-5.1.r11289-pybeef"])
for root, dirs, files in os.walk(".", topdown=False):
    for name in files:
        if name == "opt.json" or name == "opt.pckl":
            print(os.path.join(root, name))
            os.path.remove(os.path.join(root, name))
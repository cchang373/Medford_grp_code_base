# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 15:36:45 2017

@author: benjamin
"""

import sys
sys.path.insert(0,"/home/benjamin/catmap")
sys.path.insert(0,"/home/benjamin/espresso/espresso_rutile/tools")
import better_pickle_reader
import numpy as np
#plt = catmap.plt


energy_dict = better_pickle_reader.pickle_pull()

for item in energy_dict:
    if type(energy_dict[item])==list:
        print(item)
        print(energy_dict[item][0]-np.average(energy_dict[item][1]))
        #print(np.std(energy_dict[item][1]))
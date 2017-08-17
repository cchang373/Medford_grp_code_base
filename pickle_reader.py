# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 10:15:28 2017

@author: benjamin
"""

import os
from pickle import load
import numpy as np

def BEEF_rxn_eng(input_list, energy_dict):
    """
    input_list type is list, energy_dict type is dictionary. list should have two items. the fir
    st beingthe surface+adsorbate the second being the molecule alone. 4 layer rutile is assumed
    
    """
    reference = energy_dict['O2_molecules']+energy_dict['N2_molecules']+energy_dict['4_no_ads'] 
    try:
        print(input_list[1][:-10])
        print(energy_dict[input_list[0]][0]-energy_dict[input_list[1]][0]-energy_dict['4_no_ads'][0],'eV, energy from BEEF\n')
        print(np.mean(energy_dict[input_list[0]][1]-energy_dict[input_list[1]][1]-energy_dict['4_no_ads'][1]), 'eV, average energy from BEEF\n')
        print(np.std(energy_dict[input_list[0]][1]-energy_dict[input_list[1]][1]-energy_dict['4_no_ads'][1]), 'eV, error bar\n')
    except:
        print('make sure both entries are ensembles')
        

def BEEF_ref_eng(input_list, energy_dict):
    """
    input_list type is list, energy_dict type is dictionary. list should have two items. the fir
    st beingthe surface+adsorbate the second being the molecule alone. 4 layer rutile is assumed
    
    """
    reference = energy_dict['O2_molecules']+energy_dict['N2_molecules']+energy_dict['4_no_ads'] #N2 + O2 + *
    try:
        print(input_list[1][:-10])
        abs_eng = 0
        mean_eng = 0
        std_eng = 0
        ens_sum = np.array[0]*2000
        for i in len(input_list):
            abs_eng = abs_eng + energy_dict[input_list[i]][0] #abs energy eV
            ens_sum = abs_eng + energy_dict[input_list[i]][1]#eV, average energy from BEEF
        abs_eng = abs_eng - reference[0]
        ens_sum = ens_sum - reference[1]
        
    except:
        print('make sure both entries are ensembles')





#pulls in all the pickles in the pickle_barrel folder into energy_dict
#energy_dict should have the structure: {'system name': [energy from BEEF, [full 2000 member ensemble]]}
#this is also setup to handle energy dictionaries with no ensemble {'system name': energy from BEEF}
pickle_location = '/nv/hp13/bcomer3/shared/espresso_rutile/espresso_Rutile/pickle_barrel'
#pickle_location = '/home/benjamin/espresso/espresso_rutile/espresso_Rutile/pickle_barrel'
os.chdir(pickle_location)
energy_dict = {'':[0,[0]*2000]}
for root, dirs, files in os.walk(".", topdown=False):
    for name in files:
        energy_dict[name]=load(open(os.path.join(root,name),'r'))

#print(energy_dict)
#print(energy_dict['test_NO'])
BEEF_rxn_eng(['test_NO','NO_molecules'],energy_dict)
BEEF_rxn_eng(['test_N2','N2_molecules'],energy_dict)
BEEF_rxn_eng(['test_NO','NO_molecules'],energy_dict)
BEEF_rxn_eng(['test_NO2','NO2_molecules'],energy_dict)
BEEF_rxn_eng(['test_N2O','N2O_molecules'],energy_dict)
BEEF_rxn_eng(['angle_try_NO','NO_molecules'],energy_dict)

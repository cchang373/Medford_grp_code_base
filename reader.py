# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 10:15:28 2017

@author: benjamin
"""

import os
from pickle import load
import numpy as np
import sys
from ase.atoms import string2symbols
sys.path.insert(0,"/home/benjamin/catmap")
from ThemoChem import vibration_energy_ig_list, vibration_energy_list,ZPE_only
from Beef_data import beef as bf

def BEEF_ref_eng(input_list, energy_dict,multiples=None):
    """
    input_list type is list, energy_dict type is dictionary. list should be the
    names of the pickles associated with the . 
    the last item in the list surface+adsorbate. 4 layer rutile is assumed
    
    """
    if multiples==None:
        multiples = [1]*len(input_list)
    reference = []
    for i in range(2):
        reference.append(energy_dict['O2-molecules-DFT_references'][i]+energy_dict['N2-molecules-DFT_references'][i]+energy_dict['slab-slab-4_layer'][i]*multiples[-1])#N2 + O2 + *
#    try:
#        print(input_list[1][:-10])
    abs_eng = 0
    ens_sum = np.zeros((1,2000))
    for i in range(len(input_list)):
        abs_eng = abs_eng + energy_dict[input_list[i]][0]*float(multiples[i]) #abs energy eV
        ens_sum = ens_sum + energy_dict[input_list[i]][1]*float(multiples[i])#eV, average energy from BEEF
    abs_eng = abs_eng - reference[0]
    ens_sum = ens_sum - reference[1]
    std_eng = np.std(ens_sum)
    mean_eng = np.mean(ens_sum)
    return bf([abs_eng, ens_sum])
#    except:
#        print('make sure both entries are ensembles')

def get_formation_energies(energy_dict,ref_dict):
    formation_energies = {}
    for key in energy_dict: #iterate through keys
        if type(energy_dict[key])==list or isinstance(energy_dict[key]):
            E0 = [energy_dict[key][0]]
            E0.append(energy_dict[key][1])#raw energy for ensemble    
            name,details,site = key.split('-') #split key into name/site
            if 'slab' not in name: #do not include empty site energy (0)
                if site == '4_layer' or site == '2_layer_runs':
                    E0[0] -= ref_dict[site][0]
                    E0[1] -= ref_dict[site][1]#subtract slab energy if adsorbed
                #remove - from transition-states
                formula = name
                #get the composition as a list of atomic species
                try:
                    composition = string2symbols(formula)
                 #for each atomic species, subtract off the reference energy
                    for atom in composition:
                        E0[0] = E0[0] -ref_dict[atom][0]
                        E0[1] = E0[1] -ref_dict[atom][1]
                    #round to 3 decimals since this is the accuracy of DFT
        #            E0 = [np.round(x,3) for x in E0] #tried to get it to round each item in list
                    formation_energies[key] = E0
                except:
                    pass
                del E0
    return formation_energies




#pulls in all the pickles in the pickle_barrel folder into energy_dict
#energy_dict should have the structure: {'system name': [energy from BEEF, [full 2000 member ensemble]]}
#this is also setup to handle energy dictionaries with no ensemble {'system name': energy from BEEF}
def pickle_pull():
    pickle_location = '/nv/hp13/bcomer3/shared/espresso_rutile/espresso_Rutile/pickle_barrel'
    pickle_location_local = '/home/benjamin/espresso/espresso_rutile/espresso_Rutile/pickle_barrel'
    cur_dir = os.path.abspath('.')
    try:
        os.chdir(pickle_location)
    except:
            try:
                os.chdir(pickle_location_local)
            except:
                print('couldn\'t get to pickle barrel folder')
    energy_dict = {}
    for root, dirs, files in os.walk(".", topdown=False): 
        for name in files:
            energy_dict[name]=bf(load(open(os.path.join(root,name),'r')))
    os.chdir(cur_dir)
    return(energy_dict)

def vib_pickle_pull():
    pickle_location = '/nv/hp13/bcomer3/shared/espresso_rutile/espresso_Rutile/vib_pickle_barrel'
    pickle_location_local = '/home/benjamin/espresso/espresso_rutile/espresso_Rutile/vib_pickle_barrel'
    cur_dir = os.path.abspath('.')
    try:
        os.chdir(pickle_location)
    except:
            try:
                os.chdir(pickle_location_local)
            except:
                print('couldn\'t get to pickle barrel folder')
    vib_dict = {}
    for root, dirs, files in os.walk(".", topdown=False): 
        for name in files:
            vib_dict[name]=load(open(os.path.join(root,name),'r'))
    os.chdir(cur_dir)
    return(vib_dict)
    
def add_vib_energy(energy_dict,Temp,Press,spec_case=[],spec_press=[]):
    vib_dict = vib_pickle_pull()
    vib_eng_dict = {}
    has_vib = []
    #print vib_dict
    for key in vib_dict:
        if key in spec_case:
            numb = spec_case.index(key)
            vib,name,details,site = key.split('-')
            vib_eng_dict[key]=vibration_energy_ig_list(name,vib_dict[key],Temp,spec_press[numb])
        elif 'DFT' in key:
            vib,name,details,site = key.split('-')
            vib_eng_dict[key]=vibration_energy_ig_list(name,vib_dict[key],Temp,Press)
        else:
            vib_eng_dict[key]=vibration_energy_list(vib_dict[key],Temp)
        _,wo_vib = key.split('-',1)
        if wo_vib in energy_dict.keys(): #and type(energy_dict[wo_vib])==list:
            energy_dict[wo_vib][0] = energy_dict[wo_vib][0] + vib_eng_dict[key]
            energy_dict[wo_vib][1] = energy_dict[wo_vib][1] + vib_eng_dict[key]
            has_vib.append(wo_vib)
            
    for key in energy_dict.keys():
        if (key not in has_vib) and (key!='slab-slab-4_layer') and (key!='slab-slab-d_4_layer'):
            del energy_dict[key]
    return energy_dict
    
def add_ZPE(energy_dict):
    vib_dict = vib_pickle_pull()
    #print(vib_dict.keys())
    vib_eng_dict = {}
    #print vib_dict
    has_vib = []
    for key in vib_dict:
        #print(key)
        if 'DFT' in key:
            vib,name,details,site = key.split('-')
            vib_eng_dict[key]=ZPE_only(vib_dict[key])
        else:
            vib,name,details,site = key.split('-')
            vib_eng_dict[key]=ZPE_only(vib_dict[key])
        _,wo_vib = key.split('-',1)
        #print(vib_eng_dict[key])
        if wo_vib in energy_dict.keys() and type(energy_dict[wo_vib])==list:
            energy_dict[wo_vib][0] = energy_dict[wo_vib][0] + vib_eng_dict[key]
            energy_dict[wo_vib][1] = energy_dict[wo_vib][1] + vib_eng_dict[key]
            has_vib.append(wo_vib)
            #print(key)
    for key in energy_dict.keys():
        if (key not in has_vib) and (key!='slab-slab-4_layer')and (key!='slab-slab-d_4_layer'):
            del energy_dict[key]
    return energy_dict


energy_dict = pickle_pull()
add_vib_energy(energy_dict,300,101325)

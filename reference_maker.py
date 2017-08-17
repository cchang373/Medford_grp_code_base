# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 15:26:39 2017

@author: benjamin
"""
import numpy as np
from Beef_data import beef as bf
from Pourbaix import oxidation_state
from ase.atoms import string2symbols

def convert_formation_energies(energy_dict,atomic_references,composition_dict,potential):
    """
    Convert dictionary of energies, atomic references and compositions into a dictionary of formation energies
    :param energy_dict: Dictionary of energies for all species.
                        Keys should be species names and values
                        should be energies.
                        
    :type energy_dict: dict
    :param atomic_references: Dictionary of atomic reference compositions (e.g. {H2O:{H:2,O:2}})
    :type atomic_references: dict
    :param composition_dict: Dictionary of compositions
    :type composition_dict: dict
    .. todo:: Explain the keys and values for energy_dict, atomic_references, and composition_dict
    """
#    ele =-1.602176620898*10**-19
    n = len(atomic_references)
    R = np.zeros((n,n))
    e = []
    ref_offsets = {}
    atoms = sorted(atomic_references)
#    energy_dict['H2'][0] = energy_dict['H2'][0]+2*potential
#    energy_dict['H2'][1] = energy_dict['H2'][1]+2*potential
    for i,a in enumerate(atoms):
        composition = composition_dict[atomic_references[a]]
        #print(type(composition))
        e.append(energy_dict[atomic_references[a]])
        for j,a in enumerate(atoms):
            n_a = composition.get(a,0)
            R[i,j] = n_a
    if not np.prod([R[i,i] for i in range(0,n)]):
        raise ValueError('Reference set is not valid.')
    e1 = []
    e2 = []
    for i in range(len(e)):
        e1.append(e[i][0])
        e2.append(e[i][1])
    e1 = np.array(e1)
    e2 = np.array(e2)
    try:
        R_inv = np.linalg.solve(R,np.eye(n))
    except np.linalg.linalg.LinAlgError:
        raise ValueError('Reference set is not valid.')
    x = list(np.dot(R_inv,e1))
    y = list(np.dot(R_inv,e2))
    for i,a in enumerate(atoms):
        ref_offsets[a] = [x[i],y[i]]
    for key in ref_offsets:
        ref_offsets[key] = bf(ref_offsets[key])
#    new_data = {}
#    for key in energy_dict:
#        composition = composition_dict[key]
#        E1 = energy_dict[key][0]
#        E2 = energy_dict[key][1]
#        for symb in composition:
#            E1 -= ref_offsets[symb][0]*composition[symb]
#            E2 -= ref_offsets[symb][1]*composition[symb]
#        new_data[key] = [E1,E2]
        
    return ref_offsets
def rxn_energy(rxn_states, initial_state, energy_dict):
    rxn_states = [s for s in rxn_states] #make a copy
    initial_state = [[s for s in initial_state]] #make a copy
    E_p1 = []
    E_p2 = []
    for state in initial_state+rxn_states:
        for i, sp in enumerate(state): #pull off the multiplication
            if '*' in sp:
                n,s = sp.split('*')
                n = float(n)
            else:
                n = 1
                s = sp
            state[i] = [n,s]
    for i,state in enumerate(initial_state+rxn_states):
        E_p1.append(sum([d*energy_dict[p][0] for d,p in (initial_state+rxn_states)[i]]))
        E_p2.append(sum([d*energy_dict[p][1] for d,p in (initial_state+rxn_states)[i]]))
    E_r1 = sum([s*energy_dict[r][0] for s,r in initial_state[0]])
    E_r2 = sum([s*energy_dict[r][0] for s,r in initial_state[0]])
    dE1 = E_p1-E_r1
    dE2 = E_p2-E_r2
    state_energy = [list(a) for a in zip(dE1,dE2)]
    return state_energy

def get_formation_energies(energy_dict1,ref_dict,pot):
    pot= -pot#-1.05727231764
    formation_energies = {}
    for key in energy_dict1: #iterate through keys
        energy_dict = dict(energy_dict1)
        if type(energy_dict[key])==list or isinstance(energy_dict[key],bf):
            E0 = [energy_dict[key][0]]
            E0.append(energy_dict[key][1].copy())#raw energy for ensemble    
            name,details,site = key.split('-') #split key into name/site
#            if 'slab' not in name: #do not include empty site energy (0)
            if site == '4_layer' or site =='d_4_layer': #or site == '2_layer_runs':
                E0[0] -= energy_dict['slab-slab-'+site][0]
                E0[1] -= energy_dict['slab-slab-'+site][1]#subtract slab energy if adsorbed    
            formula = name
            #get the composition as a list of atomic species
            try:
                composition = string2symbols(formula)
                oxy = oxidation_state(formula)
                #print(key)
                E0[0] = E0[0] + oxy*pot
                E0[1] = E0[1] + oxy*pot
             #for each atomic species, subtract off the reference energy
                for atom in composition:
                    E0[0] = E0[0] -ref_dict[atom][0]
                    E0[1] = E0[1] -ref_dict[atom][1]
                    if site == '4_layer' and atom =='H':
                        E0[0] = E0[0] 
                        E0[1] = E0[1] 
                formation_energies[key] = bf(E0) #+oxy*pot
            except:
                if 'slab' in name:
                   formation_energies[key] = bf(E0)
            del E0
    #formation_energies['e'] = [-pot,np.zeros(2000)-pot]
    fe = formation_energies.copy()
#    SHE = 0.5*fe['H2-molecules-DFT_references'][0]-pot
#    print('H2 evolution: ')
#    print(0.5*fe['H2-molecules-DFT_references'][0]-pot)
#    print('O2 evolution: ')
#    print(-fe['H2O-molecules-DFT_references'][0]+0.5*fe['O2-molecules-DFT_references'][0]+fe['H2-molecules-DFT_references'][0]+2*pot)
#    print('NO evolution: ')
#    print(2*fe['NO-molecules-DFT_references'][0]-2*fe['H2O-molecules-DFT_references'][0]-fe['N2-molecules-DFT_references'][0])
    return formation_energies
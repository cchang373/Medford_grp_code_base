# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 14:30:48 2017

@author: benjamin
"""
import numpy as np
from glob import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os


def get_BEEF_ensemble(log):
    f = open(log)
    txt = f.read()
    f.close()
    _,E_total = txt.rsplit('total energy              =',1)
    E_total,_ = E_total.split('Ry',1)
    E_total = float(E_total.strip())
    E_total *= 13.605698
    try:
        _, ens = txt.rsplit('BEEFens 2000 ensemble energies',1)
        ens,_ = ens.split('BEEF-vdW xc energy contributions',1)
        ens.strip()
        ens_ryd = []
        for Ei in ens.split('\n'):
            if Ei.strip():
                Ei = float(Ei.strip())
        ens_ryd.append(Ei)
        ens_ryd = np.array(ens_ryd)
        ens_eV = ens_ryd * 13.605698
    except:
        ens_eV = np.array(['a'])

    return E_total, np.array(E_total + ens_eV)

def rxn_energy(products, reactants, energy_dict):
    products = [s for s in products] #make a copy
    reactants = [s for s in reactants] #make a copy
    for state in [products,reactants]:
        for i, sp in enumerate(state):
            if '*' in sp:
                n,s = sp.split('*')
                n = float(n)
            else:
                n = 1
                s = sp
            state[i] = [n,s]
    E_p = sum([s*energy_dict[p] for s,p in products])
    E_r = sum([s*energy_dict[r] for s,r in reactants])
    dE = E_p - E_r
    return dE

def get_error(products, reactants, calc_dict, true_dict):
    dE_calc = rxn_energy(products, reactants, calc_dict)
    dE_true = rxn_energy(products, reactants, true_dict)
    return np.abs(dE_calc - dE_true)

ensembles = {}
for root, dirs, files in os.walk(".", topdown=False):
    for name in files:
        if name == "log":
                ensembles[name] = get_BEEF_ensemble(os.sep())[-1]

            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            

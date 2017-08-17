# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 12:17:55 2017

@author: benjamin

"""

from ase.thermochemistry import HarmonicThermo, IdealGasThermo
import os
from ase.build import molecule
import numpy as np
#from scipy import optimize
#import numpy as np
#import matplotlib.pyplot as plt

form_dict = {
             'N2':{'atoms':molecule('N2'),'geometry':'linear','natoms':2,'symmetrynumber':2,'spin':0},
             'NO':{'atoms':molecule('NO'),'geometry':'linear','natoms':2,'symmetrynumber':1,'spin':0.5},
             'O2':{'atoms':molecule('O2'),'geometry':'linear','natoms':2,'symmetrynumber':2,'spin':1},
             'NH3':{'atoms':molecule('NH3'),'geometry':'nonlinear','natoms':4,'symmetrynumber':9,'spin':0},
             'H2':{'atoms':molecule('H2'),'geometry':'linear','natoms':2,'symmetrynumber':2,'spin':0},
             'H2O':{'atoms':molecule('H2O'),'geometry':'nonlinear','natoms':3,'symmetrynumber':2,'spin':0},
             'NO2':{'atoms':molecule('NO2'),'geometry':'nonlinear','natoms':3,'symmetrynumber':2,'spin':0},
             'N2O':{'atoms':molecule('N2O'),'geometry':'nonlinear','natoms':3,'symmetrynumber':2,'spin':0}
                }
#references at the bottom
antiones_dict={'H2O':{'A':6.20963,'B':2354.731	,'C':7.559}   #293K - 343K Gubkov, Fermor, et al., 1964 (1)

}

def vap_press(compound,Temp):
    """
    takes Temps in Kelvin
    returns Pressures in bar
    VERY APPROXIMATE
    enter compound name as a string, it should appear as it does in the 
    antiones_dict variable above. If it isn't in the dictionary feel free to 
    add it. The data is likely available in the NIST web book. Make sure you
    include the citation.
    """
    if compound not in antiones_dict.keys():
        print("vapor pressure data not available for this compound")
    else:
        A = antiones_dict[compound]['A']
        B = antiones_dict[compound]['B']
        C = antiones_dict[compound]['C']
        P = 10**(A-(B/(Temp+C)))
        return P
    
    
    
def vib_parse(path):
    """
    parses vibrational summary files from ASE
    """
    f = open(path)
    raw = f.read().splitlines()
    f.close()
    vib_eng = [] #meV out of the file, converted to eV below
    vib_freq = [] #cm^-1
    kT = 8.6173303*10**-5*300 #eV
    for line in raw:
       if '-' not in line and 'eV' not in line:
           d = line.strip()
           d,freq = d.rsplit(' ',1)
           try:
               vib_freq.append(float(freq))
               d = d.strip()
               d,eng = d.rsplit(' ',1)
               vib_eng.append(float(eng)/10**3) #converting to eV
           except:
               vib_eng.append('insert_3kT')
               vib_freq.append('insert_3kT')
    return vib_eng,vib_freq
    
    
def vibration_energy(path,temp):
    vib_eng,vib_freq = vib_parse(path)
    vib_eng = [x for x in vib_eng]
    num_3kT = len([x for x in vib_eng if x=='insert_3kT'])
    vib_eng = [x for x in vib_eng if type(x)==float]
    vib_eng.sort()
    vib_thermo = HarmonicThermo(vib_eng)
    free_eng = vib_thermo.get_helmholtz_energy(temp,verbose=False)
    free_eng = free_eng + 8.6173303*10**-5*temp*num_3kT
    return free_eng
    
def ZPE_only(vib_in):
    if not type(vib_in)==list:
        return 0
    [vibrations_list,frequency] = vib_in
    vib_eng = [x for x in vibrations_list]
    vib_eng.sort()
    for i,item in enumerate(vib_eng):
        if item =='insert_3kT':
            vib_eng[i] = 30/8065.54429
    vib_eng = [x for x in vib_eng if type(x)==float]
    vib_thermo = HarmonicThermo(vib_energies=vib_eng)
    free_eng = vib_thermo.get_ZPE_correction()
    #print vib_thermo.get_gibbs_energy(temp,pressure,verbose=True)
    return free_eng
    
def vibration_energy_ig_list(formula ,vib_in,temp,pressure):
    if not type(vib_in)==list:
        return 0
    [vibrations_list,frequency] = vib_in
    vib_eng = [x for x in vibrations_list]
    vib_eng.sort()
    for i,item in enumerate(vib_eng):
        if item =='insert_3kT':
            vib_eng[i] = 30/8065.54429
    vib_eng = [x for x in vib_eng if type(x)==float]
    if formula in form_dict.keys():
        vib_thermo = IdealGasThermo(vib_energies=vib_eng,**form_dict[formula])
        free_eng = vib_thermo.get_gibbs_energy(temp,pressure,verbose=False)
        #print vib_thermo.get_gibbs_energy(temp,pressure,verbose=True)
        return free_eng
    else:
        return 0

def vibration_energy_list(vib_in,temp):
    if not type(vib_in)==list:
        return 0
    [vibrations_list,frequency] = vib_in
    vib_eng = [x for x in vibrations_list]
    for i,item in enumerate(vib_eng):
        if item =='insert_3kT':
            vib_eng[i] = 30/8065.54429
    vib_eng = [x for x in vib_eng if type(x)==float]
    vib_thermo = HarmonicThermo(vib_eng)
    free_eng = vib_thermo.get_helmholtz_energy(temp,verbose=False)
    return free_eng
    

        
def gas_vibrations(path,temp,pressure):
    path2 = path
    path_list = path2.split('/')
    formula = path_list[-3]
    if formula in form_dict.keys():
        vib_eng,vib_freq = vib_parse(path)
        vib_eng = [x for x in vib_eng if type(x)==float]
        #vib_eng = np.array(vib_eng)
        #vib_eng = [np.max(vib_eng)]
        #print vib_thermo.get_gibbs_energy(temp,pressure,verbose=True)

        return [vib_eng,vib_freq]
    else:
        return 0

def surface_vibrations(path,temp):
    vib_eng,vib_freq = vib_parse(path)
    return [vib_eng,vib_freq]

#I never could get these to work, work on these in the future
def find_3kT_eng(temp):
    eng = optimize.bisect(vib_eng_3kT,0.00001,0.5,args=(temp))
    return eng

def vib_eng_3kT(energy,T):
    k = 8.6173303*10**-5
#    print np.exp(energy/kT)
    err = energy/(np.exp(energy/k/T)-1)-3*k*T
    return err


def Gibbs_isotherm(energy_dict,index,temp=300,area=6.5797272*5.96159647):
    """
    takes in a dictionary of energies in the form [energy, ensemble] and 
    returns the coverage from a gibs isotherm at the temp and pressure
    """
    prob = {}
    k = 8.617e-5
    T = temp
    a = area
    for key in energy_dict.keys():
        new_entry = energy_dict[1][key][index]
        new_entry = np.append(new_entry,state_free_eng_dict[key][index])
        state_free_eng_dict[key] = np.mean(new_entry)
    for item in energy_dict.copy():
        prob[item]=np.exp(-state_free_eng_dict[item].copy()/(k*T))
    tot = 0
    for item in state_free_eng_dict.copy():
            tot += prob[item]
    tot_P = tot
    check = 0
    for item in state_free_eng_dict.copy():   
        prob[item] = prob[item]/tot_P
        check += prob[item]

"""
References:
1. Gubkov, A.N.; Fermor, N.A.; Smirnov, N.I., Vapor Pressure of Mono-Poly 
        Systems, Zh. Prikl. Khim. (Leningrad).


"""
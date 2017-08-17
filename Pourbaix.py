# -*- coding: utf-8 -*-
"""
Created on Thu May  4 16:48:38 2017

@author: benjamin
"""
from ase.atoms import string2symbols
import numpy as np
import matplotlib.cm as cm
import pylab as plt
from Beef_data import beef as bf

k = 8.6173303e-5
T = 300

def N_oxidation_state(compound):
    """
    finds the oxidation state of compounds with the formual NxHyOz. 
    It's not remotely robust
    """
    cmpd = string2symbols(compound)
    numN = len([x for x in cmpd if x=='N'])
    os = 0
    for atom in cmpd:
        if atom == 'O':
            os = os-2
        elif atom == 'H':
            os = os+1
    N_os = -os/numN
    return N_os

def oxidation_state(compound):
    """
    finds the oxidation state of compounds with the formual NxHyOz. 
    It's not remotely robust
    """

    if compound=='slab': #or compound == 'O2' or compound == 'H2' or compound=='O4':
        return 0
    os = 0
    cmpd = string2symbols(compound)
    for atom in cmpd:
        if atom == 'O':
            os = os-2
        elif atom == 'H':
            os = os+1
    os = -os
    return os

def num_hydrogens(compound):
    try:
        cmpd = string2symbols(compound)
        numH = len([x for x in cmpd if x=='H'])
    except:
        numH = 0
    return numH
def num_nitrogens(compound):
    try:
        cmpd = string2symbols(compound)
        numN = len([x for x in cmpd if x=='N'])
    except:
        numN = 0
    return numN

def pH_correct(oxy,numH,pH,Temp):
    eng = k*np.log(10)*pH*oxy*Temp
    return eng

def potential_correct(os,potential):
    eng = potential*os
    return eng

def Pourbaix(energy_dict,pH_range,Potential,Temp):
    """
    Returns a matrix with constant pH along the rows and constant 
    """
    n = 0
    surf_eng = {}
    for key in energy_dict: #iterate through keys
        other = '2_layer_runs'not in key and 'rect' not in key and key!='e'
        if (type(energy_dict[key])==list or isinstance(energy_dict[key],bf)) and other:
            energy_dict = dict(energy_dict)
            E0 = energy_dict[key][0]
    #        E0.append(energy_dict[key][1].copy())#raw energy for ensemble    
            name,details,site = key.split('-') #split key into name/site
            if site =='4_layer':
                surf_eng[key]=np.zeros((len(pH_range.tolist()),len(Potential.tolist())))
                #print(surf_eng[key].size)
                for i,pH in enumerate(pH_range):
                    for j,U in enumerate(Potential.tolist()):
                        numH = num_hydrogens(name)
                        oxy = oxidation_state(name)
                        E0 = E0 - potential_correct(oxy,U) - pH_correct(oxy,numH,pH,Temp)
        #                E0[1] = (E0[1]-potential_correct(0,oxy,U) +pH_correct(numH,pH))
                        surf_eng[key][i,j]=E0
                        n = n+1
                        #print(n)
                        E0 = energy_dict[key][0]
    return surf_eng

def Pourbaix_test(energy_dict,Potential,pH_range,Temp):
    surf_eng = {}
    n=0
    for key in energy_dict: #iterate through keys
        energy_dict = dict(energy_dict)
        E0 = energy_dict[key][0]
        name,details,site = key.split('-') #split key into name/site
#        if 'slab' not in name: #do not include empty site energy (0)
        surf_eng[key]=np.zeros((len(pH_range.tolist()),len(Potential.tolist())))

        for i,pH in enumerate(pH_range):
            for j,U in enumerate(Potential.tolist()):
    #            oxy = N_oxidation_state(name)
                numH = num_hydrogens(name)
                oxy = oxidation_state(name)
                E0 = (E0-(potential_correct(oxy,U) +pH_correct(oxy,numH,pH,Temp))*energy_dict[key][2])/energy_dict[key][1]
                #E0 = (E0)/energy_dict[key][1]-potential_correct(oxy,U) -pH_correct(oxy,numH,pH,Temp)
                #print(n)
                n = n+1
                surf_eng[key][i,j]=E0
                E0 = energy_dict[key][0] 
    return surf_eng


def Pourbaix_Diagram(surf_eng,pH_range,U_range):
    species = [x for x in surf_eng.keys()]
    pourbaix = np.empty((len(U_range),len(pH_range)),int)
    for j,U in enumerate(U_range):#find the lowest energy species at a given (u,pH)
        for i,pH in enumerate(pH_range):
            engs = 10000
            for item in surf_eng.keys():
                eng = surf_eng[item][i,j]
                if eng < engs:
                    engs = eng
                    pourbaix[j,i] = species.index(item)
    return pourbaix
    
energy_dict ={
                'Ag-metal-':[0,4,1],
                'O-1/4-ML':[2.15,4,1],
                'O-1/3-ML':[4.46,6,2],
#                'O-1/2-ML':[5.39,4,2],
#                'OH-1/2-ML':[3.28,6,3],
                'OH-1/6-ML':[0.93,6,1],
                'OH-1/3-ML':[2.03,6,2],
                #'Ag-diss-ML':[4.8,6,6]
}

"""energy_dict ={
                'Ag-P-ML':[0,4],
                'O-1/4-ML':[2.10,4],
                'O-1/3-ML':[4.36,6],
#                'O-1/2-ML':[5.39,4],
                'OH-1/2-ML':[2.23,6],
                'OH-1/6-ML':[0.58,6],
                'OH-1/3-ML':[1.33,6],
}"""


"""
pH_range = np.linspace(0,14,20)
U_range = np.linspace(-1,1.5,20)
surf_eng = Pourbaix_test(energy_dict,U_range,pH_range,300)
pourbaix = Pourbaix_Diagram(surf_eng,pH_range,U_range)
#np.set_printoptions(threshold='nan')
#print pourbaix[::-1]
#np.set_printoptions(threshold='nan')
#print(pourbaix)
#ax = plt.gca()
#ax.imshow(pourbaix, cmap=cm.Accent,
#          extent=[min(pH_range), max(pH_range),min(U_range), max(U_range)],
#          origin='lower',
#          aspect='auto')
#ax.legend()
#print surf_eng
for item in surf_eng:
    plt.plot(U_range,surf_eng[item][0,:],label=item)
plt.legend()
#plt.show()
plt.savefig('pourbaix2.jpg')"""
#

# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 10:01:38 2017

@author: benjamin
"""
import sys
import reader
import reference_maker
from reference_maker import get_formation_energies
from plotter import move_free_energy,plot_FED,rxn_energy,eng_add_OH_hack,plot_FED_return
from Pourbaix import *
#import catmap.data.parameter_data#.experimental_gas_formation_energies as egfe
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from ThemoChem import vap_press
import pickle

fig_dump = '/home/benjamin/paper1'

P_water = vap_press('H2O',300) #bar
P_N2 = 0.8 #atm
P_O2 = 0.2 #atm

energy_dict = reader.pickle_pull()
#print energy_dict['NO-molecules-DFT_references'][0]
#print(2*energy_dict['NO-molecules-DFT_references'][0]-energy_dict['O2-molecules-DFT_references'][0]-energy_dict['N2-molecules-DFT_references'][0])

energy_dict = reader.add_vib_energy(energy_dict,300,101325,
                                                  spec_case=['vib-H2O-molecules-DFT_references',
                                                             'vib-N2-molecules-DFT_references',
                                                             'vib-O2-molecules-DFT_references',
                                                             
                                                             ],
                                                  spec_press=[P_water*100000,
                                                              P_N2*101325,
                                                              P_O2*101325
                                                  
                                                  ]
)
#print(np.mean(energy_dict['slab-slab-4_layer'][1]-energy_dict['slab-slab-d_4_layer'][1]-0.5*energy_dict['O2-molecules-DFT_references'][1]))
energy_dict = eng_add_OH_hack(energy_dict)
atomic_references = {'N':'N2','O':'O2','H':'H2O','C':'CO2'}
composition_dict = {'N2':{'N':2},'NO':{'N':1,'O':1},'O2':{'O':2},
                    'N2O':{'N':2,'O':1},'N2O2':{'N':1,'O':2},
                    'NO2':{'N':1,'O':2},'CO2':{'C':1,'O':2},
                    'OH':{'H':1,'O':1},
                    'H2':{'H':2}, 'H2O':{'H':2,'O':1},'NH3':{'N':1,'H':3},
                    'CO':{'C':1,'O':1},
                    'CH4':{'C':1,'H':4}
                    }
compound_en_dict = {}
for key in energy_dict.copy():
    name,details,site = key.split('-') #split key into name/site
    if site == 'DFT_references' and 'molecules_2' not in key:#and (('N' in name or 'O' in name) and 'H' not in name and 'C' not in name):
        compound_en_dict[name] = [energy_dict[key][0]]
        compound_en_dict[name].append(energy_dict[key][1].copy())
      
#pot = 0.00833
pot = 0
surface_area = 6.5797272*5.96159647
ref_offsets = reference_maker.convert_formation_energies(compound_en_dict.copy(),atomic_references.copy(),composition_dict.copy(),pot)
#ref_offsets = ref_to_OH_hack(ref_offsets.copy(),energy_dict.copy())
formation_energies = get_formation_energies(energy_dict.copy(),ref_offsets.copy(),pot)
#print(2*energy_dict['NO-molecules-DFT_references'][0]+2*energy_dict['H2-molecules-DFT_references'][0]-2*energy_dict['H2O-molecules-DFT_references'][0]-energy_dict['N2-molecules-DFT_references'][0])
#print(-energy_dict['H2O-molecules-DFT_references'][0]+0.5*energy_dict['O2-molecules-DFT_references'][0]+energy_dict['H2-molecules-DFT_references'][0])
formation_energies= eng_add_OH_hack(formation_energies)
fe = formation_energies.copy()

#print(fe['H2-molecules-DFT_references'][0]+2*fe['OH-H2_H2O_hacked-DFT_references'][0]-2*fe['H2O-molecules-DFT_references'][0])
#print(fe['H2O-molecules-DFT_references'][0]-fe['H2-molecules-DFT_references'][0]-0.5*fe['O2-molecules-DFT_references'][0])
#print(fe['H2O-molecules-DFT_references'][0]-0.5*fe['O2-molecules-DFT_references'][0]-fe['H2-molecules-DFT_references'][0])
#print('H2 evolution: ')
#print(0.5*fe['H2-molecules-DFT_references'][0]+pot)
#print('O2 evolution: ')
#print((-fe['H2O-molecules-DFT_references'][0]+0.5*fe['O2-molecules-DFT_references'][0]-2*pot)/2)

fig, ax = plt.subplots(1,1)
initial_state = ['slab-slab-4_layer','N2-molecules-DFT_references']

states = [
            ['N2-test-4_layer'],
            ['2*N-over_O-4_layer'],
#            ['2*NH-test-4_layer'],
#            ['2*NH2-test-4_layer'],
#            ['2*NH3-test-4_layer'],
#            ['2*NH3-molecules-DFT_references']


]
labels = ['N2','N2*','2N*',]
s = rxn_energy(states,initial_state,formation_energies.copy())
f = open('N_diss.pckl','w+')
pickle.dump(s,f)
#print(s)
p1 = plot_FED_return(s,[-15]*2,labels,pot,atomic_references,ax)

del s
states = [
            ['N2O-test-4_layer','O-test-4_layer'],
            ['2*NO-test-4_layer'],
#            ['2*NH2-test-4_layer'],
#            ['2*NH3-test-4_layer'],
#            ['2*NH3-molecules-DFT_references']
            ]
labels = ['N2','N2O*+O*','2NO*',]

s = rxn_energy(states,initial_state,formation_energies.copy())
#print(s)
p2 = plot_FED_return(s,[-15]*2,labels,pot,atomic_references,ax)
del s
states = [
            ['H2NNH2-test-4_layer'],
            ['2*NH2-test-4_layer'],
#            ['2*NH2-test-4_layer'],
#            ['2*NH3-test-4_layer'],
#            ['2*NH3-molecules-DFT_references']
            ]
labels = ['N2','N2H4*','2NH2*',]
s = rxn_energy(states,initial_state,formation_energies.copy())
#print(s)
p3 = plot_FED_return(s,[-15]*2,labels,pot,atomic_references,ax)

#plt.savefig('zero_pot.pdf')
#fig = plt.figure()
#
#surface_area = 6.5797272*5.96159647
##ref_test, ref_offsets = reference_maker.convert_formation_energies(compound_en_dict.copy(),atomic_references,composition_dict,pot)
##formation_energies = get_formation_energies(energy_dict.copy(),ref_offsets.copy(),pot)
#initial_state = ['slab-slab-4_layer','N2-molecules-DFT_references']
#states2 = [
#            ['N2-test-4_layer',],
#            ['N2H-bridge_Ti-4_layer',],
#            ['HNNH-test-4_layer',],
#            ['H2NNH-test-4_layer',],
#            ['H2NNH2-test-4_layer',],
#            ['NH3-test-4_layer','NH2-test-4_layer',],
#            ['2*NH3-test-4_layer',],
#            ['2*NH3-molecules-DFT_references',],
#
#]
#del s
##fig = plt.figure()
#
##plt.subplots_adjust(hspace=0.001)
#labels = ['N$_2$','N$_2$*','HNN*','HNNH*','H$_2$NNH*','H$_2$NNH$_2$*','NH$_3$*+NH$_2$*','2NH$_3$*','2NH$_3$']
#s = rxn_energy(states2,initial_state,formation_energies.copy())
#
##plot_FED(s,[-15]*8,labels,pot,fig_dump+'associative_path_bg',atomic_references)
############################
#initial_state = ['slab-slab-4_layer','N2-molecules-DFT_references']
#states2 = [
#            ['N2-test-d_4_layer'],
#            ['N2H-test-d_4_layer',],
#            ['HNNH-test-d_4_layer',],
#            ['HNNH2-test-d_4_layer',],
#            ['H2NNH2-test-d_4_layer',],
#            ['NH3-test-4_layer','NH2-over_vac_2_int-d_4_layer',],
#            ['NH3-test-4_layer','NH3-test-d_4_layer',],
#            ['2*NH3-molecules-DFT_references',],
#
#]
#del s
#labels = ['N$_2$','N2$_*$','N$_2$H*','HNNH*','HNNH$_2$*','H$_2$NNH$_2$*','NH$_2$+NH$_3$*','2NH$_3$*','2NH$_3$']
#s = rxn_energy(states2,initial_state,formation_energies.copy())
#print(s)
#plot_FED(s,[-15]*8,labels,pot,fig_dump+'defect_reduction_bg',atomic_references)
#
############################
#initial_state = ['slab-slab-4_layer','N2-molecules-DFT_references']
#states2 = [
#            ['N2-test-4_layer'],
#            ['N2H-test-4_layer',],
#            ['HNNH-bridge_Ti-4_layer',],
#            ['H2NNH-test-4_layer',],
#            ['H2NNH2-test-4_layer',],
#            ['NH2-test-4_layer','NH3-test-4_layer'],
#            ['2*NH3-test-4_layer'],
#            ['2*NH3-molecules-DFT_references',],
#
#]
#del s
#labels = ['N$_2$','N2$_*$','N$_2$H*','HNNH*','HNNH$_2$*','H$_2$NNH$_2$*','2NH$_3$']
#s = rxn_energy(states2,initial_state,formation_energies.copy())
##print(s)
##plot_FED(s,[-15]*6,labels,pot,fig_dump+'associative_path_eq',atomic_references)
#
#############################"""
#initial_state = ['slab-slab-4_layer','N2-molecules-DFT_references','2*H2O-molecules-DFT_references']
#states2 = [
#            ['N2-test-4_layer','2*H2O-molecules-DFT_references'],
#            ['N2-test-4_layer','OH-test-4_layer','H2O-molecules-DFT_references'],
#            ['N2-test-4_layer','O-test-4_layer','H2O-molecules-DFT_references'],
#            ['N2O-test-4_layer','H2O-molecules-DFT_references'],
#            ['N2O-test-4_layer','OH-test-4_layer'],
#            ['N2O-test-4_layer','O-test-4_layer'],
#            ['N2O2-test-4_layer'],
#            ['2*NO-test-4_layer'],
#            ['2*NO-molecules-DFT_references'],
#
#]
#del s
#labels = ['N$_2$','N$_2$*','N$_2$*+OH*','N$_2$*+O*','N$_2$O*','N$_2$O+OH*','N$_2$O+O*','N$_2$O$_2$*','2NO*','2NO']
#s = rxn_energy(states2,initial_state,formation_energies.copy())
##print(s)
##print s
##print(s)
##plot_FED(s,[-15]*9,labels,pot,fig_dump+'oxidative_O_step_eq',atomic_references)
############################
#
#initial_state = ['slab-slab-4_layer','NO-molecules-DFT_references']
#states2 = [
#            ['NO-test-4_layer'],
#            ['ONH-test-4_layer',],
#            ['HNOH-test-4_layer',],
#            ['NH2OH-test-4_layer'],
#            ['NH2-test-4_layer','H2O-molecules-DFT_references'],
#            ['NH3-test-4_layer','H2O-molecules-DFT_references'],
#            ['NH3-molecules-DFT_references','H2O-molecules-DFT_references'],
#
#]
##del s
#labels = ['NO','NO*','HNO*','HN-OH*','H$_2$N-OH*','NH$_2$*','NH$_3$*','NH$_3$']
#s = rxn_energy(states2,initial_state,formation_energies.copy())
##print(s)
##plot_FED(s,[-15]*7,labels,pot,fig_dump+'NO_NH3_eq',atomic_references)
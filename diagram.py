# -*- coding: utf-8 -*-
"""
Created on Fri May  5 15:58:45 2017

@author: benjamin
"""
import sys
import reader
import reference_maker
from reference_maker import get_formation_energies
from plotter import move_free_energy,plot_surface_energy_diagram,Probability,Coverage,strip_high_energy,get_formation_energies_photochemistry,eng_add_OH_hack
from Pourbaix import *
#import catmap.data.parameter_data#.experimental_gas_formation_energies as egfe
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from ThemoChem import vap_press
k = 8.617330350*10**-5
T = 300
fig_dump = '/home/benjamin/paper1/'

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
atomic_references = {'N':'N2','O':'O2','H':'H2O','C':'CO2'}
composition_dict = {'N2':{'N':2},'NO':{'N':1,'O':1},'O2':{'O':2},
                    'N2O':{'N':2,'O':1},'N2O2':{'N':1,'O':2},
                    'NO2':{'N':1,'O':2},'CO2':{'C':1,'O':2},
                    'H2':{'H':2}, 'H2O':{'H':2,'O':1},'NH3':{'N':1,'H':3},
                    'CO':{'C':1,'O':1},
                    'CH4':{'C':1,'H':4}
                    }
compound_en_dict = {}
eng_add_OH_hack(energy_dict)
for key in energy_dict.copy():
    name,details,site = key.split('-') #split key into name/site
    if site == 'DFT_references' and 'molecules_2' not in key:#and (('N' in name or 'O' in name) and 'H' not in name and 'C' not in name):
        compound_en_dict[name] = [energy_dict[key][0]]
        compound_en_dict[name].append(energy_dict[key][1].copy())


pot = 0
surface_area = 6.5797272*5.96159647
band_details = [3.03,2.9]

ref_offsets = reference_maker.convert_formation_energies(compound_en_dict.copy(),atomic_references,composition_dict,pot)
#formation_energies = get_formation_energies_photochemistry(energy_dict.copy(),ref_offsets.copy(),band_details)
formation_energies = get_formation_energies(energy_dict.copy(),ref_offsets.copy(),pot)
for item in formation_energies.keys():
    if abs(formation_energies[item][0])>15:
        del energy_dict[item]
#formation_energies = get_formation_energies_photochemistry(energy_dict.copy(),ref_offsets.copy(),band_details)
formation_energies = get_formation_energies(energy_dict.copy(),ref_offsets.copy(),pot)
for item in formation_energies.keys():
    if abs(formation_energies[item][0])<50:
        print('\''+item+'\'\t'+str(formation_energies[item][0])+'\t'+str(np.std(formation_energies[item][1])))

eng,ens,STP,rng,P=move_free_energy(energy_dict.copy(),ref_offsets.copy(),surface_area,'2*H','H2O')
#eng,ens = strip_high_energy(eng,ens,1)
prob_H2O = Probability(eng.copy(),ens.copy(),rng,P,surface_area,STP,'H_2O')
prob_H2O.savefig(fig_dump+'H2O_prob.pdf')
cov_H2O = Coverage(eng.copy(),ens.copy(),rng,P,surface_area,STP,'H_2O')
cov_H2O.savefig(fig_dump+'H2O_coverage.pdf')
surf_fig_H2O = plot_surface_energy_diagram(eng,P,rng,STP,'H_2O')
del eng,ens,STP,rng,P
surf_fig_H2O.savefig(fig_dump+'surf_H2O.pdf')
eng,ens,STP,rng,P=move_free_energy(energy_dict.copy(),ref_offsets.copy(),surface_area,'2*N','N2')
#eng,ens = strip_high_energy(eng.copy(),ens,1)
#print(len(eng['N2-test-4_layer']))
#print(len(P))
prob_N = Probability(eng.copy(),ens.copy(),rng,P,surface_area,STP,'N_2')
prob_N.savefig(fig_dump+'N_prob.pdf')
cov_N = Coverage(eng.copy(),ens.copy(),rng,P,surface_area,STP,'N_2')
cov_N.savefig(fig_dump+'N_coverage.pdf')
surf_fig_N = plot_surface_energy_diagram(eng.copy(),P,rng,STP,'N_2')
#fig = matplotlib.figure.AxesStack
#fig.add(prob_N,1)
#fig.show()
#fig3, (ax1,ax2) = plt.subplots(1,2,sharey=True, sharex=True,figsize=(10,5))
#line3, = ax1.plot(prob_H2O.ax.get_data()[0], prob_H2O.ax.get_data()[1])
#fig3.savefig('surf.jpg')
#fig.add(prob_N,1)
surf_fig_N.savefig(fig_dump+'surf_N.pdf')


atomic_references = {'N':'N2','O':'O2','H':'H2','C':'CO2'}
composition_dict = {'N2':{'N':2},'NO':{'N':1,'O':1},'O2':{'O':2},
                    'N2O':{'N':2,'O':1},'N2O2':{'N':1,'O':2},
                    'NO2':{'N':1,'O':2},'CO2':{'C':1,'O':2},
                    'H2':{'H':2}, 'H2O':{'H':2,'O':1},'NH3':{'N':1,'H':3},
                    'CO':{'C':1,'O':1},
                    'CH4':{'C':1,'H':4}
                    }
ref_offsets = reference_maker.convert_formation_energies(compound_en_dict.copy(),atomic_references,composition_dict,pot)
formation_energies = get_formation_energies(energy_dict.copy(),ref_offsets.copy(),pot)

#print(formation_energies['O2-test-4_layer'])
U_range = np.arange(-6,3.2,0.2)
pH_range = np.arange(0,14,1)
surf_eng = Pourbaix(formation_energies.copy(),pH_range,U_range,T)
pourbaix = Pourbaix_Diagram(surf_eng,pH_range,U_range)
species = [x for x in surf_eng.keys()]
np.set_printoptions(threshold='nan')
print(pourbaix[::-1])
print(species[5])
print(species[68])
#print(species[14])
print(species[16])
U_range = np.arange(-0.9,3.2,0.01)
pH_range = np.arange(0,14,0.01)
surf_eng = Pourbaix(formation_energies,pH_range,U_range,T)

#print surf_eng
species = [x for x in surf_eng.keys()]


pourbaix = Pourbaix_Diagram(surf_eng,pH_range,U_range)
np.set_printoptions(threshold='nan')
#print pourbaix[::-1]
#print(species[27])

plt.cla()
plt.clf()
ax4 = plt.gca()
####
#plot bandgap
x1 = np.array([0,14])
y1 = np.array([-0.1+3.03,-0.1+3.03-14*k*T*np.log(10)])#-14*k*T*np.log(10)])
x2 = np.array([0,14])
y2 = np.array([-0.1,-0.1-14*k*T*np.log(10)])#-14*k*T*np.log(10)])
ax4.plot(x1,y1,linestyle=(0, (5, 10)),color='black')
ax4.plot(x2,y2,linestyle=(0, (5, 10)),color='black')
####
#plot RHE
x1 = np.array([0,14])
y1 = np.array([1.2,1.2-14*k*T*np.log(10)])
x2 = np.array([0,14])
y2 = np.array([0,-14*k*T*np.log(10)])
ax4.plot(x1,y1,linestyle=(0, (5, 10)),color='blue')
ax4.plot(x2,y2,linestyle=(0, (5, 10)),color='blue')
####
ax4.set_ylabel('Potential vs SHE (V)')
ax4.set_xlabel('pH')
ax4.imshow(pourbaix, cmap=cm.Accent,
          extent=[min(pH_range), max(pH_range),min(U_range), max(U_range)],
          origin='lower',
          aspect='auto')

#ax.legend()
    

"""for item in surf_eng:
    if abs(surf_eng[item][0,0])<=10:
        plt.plot(U_range,surf_eng[item],label = item)
        plt.savefig('Prob_diagram.jpg')"""
plt.savefig(fig_dump+'Pourbaix_crude.pdf')
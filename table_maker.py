# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:42:07 2017

@author: benjamin
"""

import sys
sys.path.insert(0,"/home/benjamin/catmap")
sys.path.insert(0,"/home/benjamin/espresso/espresso_rutile/tools")
import better_pickle_reader
import reference_maker
from plotter import get_formation_energies
from mapping import mapping_dict
#import catmap.data.parameter_data#.experimental_gas_formation_energies as egfe
import numpy as np
from ThemoChem import vap_press

#P_water = vap_press('H2O',300)

energy_dict = better_pickle_reader.pickle_pull()
#print energy_dict['NO-molecules-DFT_references'][0]
#print(2*energy_dict['NO-molecules-DFT_references'][0]-energy_dict['O2-molecules-DFT_references'][0]-energy_dict['N2-molecules-DFT_references'][0])

energy_dict = better_pickle_reader.add_ZPE(energy_dict
#                                                    ,300,101325,
#                                                  spec_case=['H2O-molecules-DFT_references'],
#                                                  spec_press=P_water*100000
)


atomic_references = {'N':'N2','O':'H2O','H':'H2','C':'CO2'}
composition_dict = {'N2':{'N':2},'NO':{'N':1,'O':1},'O2':{'O':2},
                    'N2O':{'N':2,'O':1},'N2O2':{'N':1,'O':2},
                    'NO2':{'N':1,'O':2},'CO2':{'C':1,'O':2},
                    'H2':{'H':2}, 'H2O':{'H':2,'O':1},'NH3':{'N':1,'H':3},
                    'CO':{'C':1,'O':1},
                    'CH4':{'C':1,'H':4}
                    }
compound_en_dict = {}

for key in energy_dict.copy():
    name,details,site = key.split('-') #split key into name/site
    if site == 'DFT_references':#and (('N' in name or 'O' in name) and 'H' not in name and 'C' not in name):
        compound_en_dict[name] = [energy_dict[key][0]]
        compound_en_dict[name].append(energy_dict[key][1].copy())
          
pot = 0

ref_test, ref_offsets = reference_maker.convert_formation_energies(compound_en_dict.copy(),atomic_references,composition_dict,pot)
formation_energies = get_formation_energies(energy_dict.copy(),ref_offsets.copy(),pot)
fec = formation_energies.copy()
for item in formation_energies.keys():
    if item in mapping_dict.keys():
        formation_energies[mapping_dict[item]]=formation_energies[item]
        del formation_energies[item]
    else:
        del formation_energies[item]


cap = 'Energies of all species used relative to a bare slab, H$_2$,N$_2$, and H$_2$O at 0K (ZPE included)'
f = open('/home/benjamin/pap1_2/8194080mfnpxctdkqmd/figures/table.tex','w+')
f.write('\\begin{center}\n\\begin{longtable}{|l|l|l|l|} \n\caption{'+cap+'}\n\label{tab:energies}\n\endfirsthead\n\endhead\n\hline\n')
f.write('compound & BEEF energy & ensemble mean & ensemble standard deviation  \\\ \hline\ \n')
fe = formation_energies
for item in fe.keys():    
    if abs(fe[item][0])>50:
        #print(item)
        del fe[item]
for item in fe:
    #print(item)    
#    name,detail,site = item.split('-')

    f.write(item+'&'+str(np.round(fe[item][0],2))+'&'+str(np.round(np.mean(fe[item][1]),2))+'&'+str(np.round(np.std(fe[item][1]),2))+'\\\  \n')
f.write(' \hline \end{longtable}\n\end{center}\n')
f.close()

cap = 'calculated equilibrium geometries'
f = open('/home/benjamin/pap1_2/8194080mfnpxctdkqmd/figures/pic-table.tex','w+')
f.write('\\begin{center}\n\\begin{longtable}{|l|l|l|l|}\caption{'+cap+'}\n\label{tab:pic}\n\endfirsthead\n\endhead\n\hline\n')
f.write('compound & Top View & Side View \\\ \hline\ \n')
fe = formation_energies
for item in fe.keys():    
    if abs(fe[item][0])>50:
        #print(item)
        del fe[item]
for item in fec.keys():
    if item in mapping_dict.keys() and item !='slab-slab-4_layer':
        #print(item)    
    #    name,detail,site = item.split('-')
    
        f.write(mapping_dict[item]+'&'+'\includegraphics[width=0.35\\textwidth, height=55mm]{figures/compounds/'+item+'-top.png}'
        +'&'+'\includegraphics[width=0.25\\textwidth, height=65mm]{figures/compounds/'+item+'-xside.png}'
        +'\\\ \hline\ \n')
f.write('\end{longtable}\n\end{center}')

vibs = better_pickle_reader.vib_pickle_pull()
#print vibs
cap = 'Vibrational frequencies of all species'
f = open('/home/benjamin/pap1_2/8194080mfnpxctdkqmd/figures/vib-table.tex','w+')
f.write('\\begin{center}\n\\begin{longtable}{|l|l|}\n\caption{'+cap+'}\n\label{tab:vibrations}\n\endfirsthead\n\endhead\n \n\hline\n')
f.write('compound & Vibrational Frequencies (cm$^{-1}$) \\\ \hline \n')
fe = formation_energies
for item in fe.keys():    
    if abs(fe[item][0])>50:
        #print(item)
        del fe[item]
for item in fec.keys():
    if item in mapping_dict.keys() and item !='slab-slab-4_layer':
        #print(item)    
    #    name,detail,site = item.split('-')
        N=1
        vib = ' \makecell{'
        for i in range(len(vibs['vib-'+item][1])):
            if vibs['vib-'+item][1][i] == 'insert_3kT':
                vib = vib+'30'+ ', '
            elif vibs['vib-'+item][1][i]< 30:
                vib = vib+'30'+ ', '
            else:
                vib = vib +str(vibs['vib-'+item][1][i])+ ', '
            if N % 3 ==0:
                vib = vib + ' \\\ '
            N = N+1
        vib = vib+'} '
        f.write(mapping_dict[item]+'&'+ vib
        +'\\\ \hline \n')
f.write('\end{longtable}\n\end{center}')
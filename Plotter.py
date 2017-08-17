# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 09:57:36 2017

@author: benjamin
"""
import sys
from ase.atoms import string2symbols
sys.path.insert(0,"/home/benjamin/catmap")
sys.path.insert(0,"/home/benjamin/espresso/espresso_rutile/tools")
import better_pickle_reader
import reference_maker
from catmap.analyze import MechanismPlot
#import catmap.data.parameter_data#.experimental_gas_formation_energies as egfe
import numpy as np
import pylab as plt
from ThemoChem import vap_press, Gibbs_isotherm
#from Probability import probabilities_ens
#plt = catmap.plt
file_name = 'NO_mechanism.txt'
"""
def make_ref_dict(energy_dict,ref_mapping):
    ref_dict = {}
    for en_key in energy_dict.keys():
        name,details,site = key.split('-')
        for key in ref_mapping.keys():
            if name == ref_mapping[key]:
                E0 = energy_dict[ref_mapping[key]] #raw energy
"""
def probabilities(energy_dict,ens_dict,index):
    k = 8.617e-5
    T = 300
    a = 6.5797272*5.96159647
    for key in energy_dict.keys():
        new_entry = ens_dict[key][index]
        new_entry = np.append(new_entry,energy_dict[key][index])
        energy_dict[key] = new_entry
    P_dict = {}
    for sp,ens in energy_dict.items():
        P_dict[sp] = []
        for i, e in enumerate(ens):
            x = np.sum(np.exp(-e*a/(k*T)))
            P_dict[sp] += [x]

    full_dict = {}
    for sp in P_dict:
        full_dict[sp] = 0
        for i in range(len(P_dict[sp])):
            full_dict[sp] += P_dict[sp][i]/sum([P_dict[sp_i][i] for sp_i in P_dict])
        full_dict[sp] /= len(P_dict[sp]) #normalize
    return full_dict

def get_formation_energies(energy_dict1,ref_dict,pot):
    formation_energies = {}
    for key in energy_dict1: #iterate through keys
        energy_dict = dict(energy_dict1)
        if type(energy_dict[key])==list:
            E0 = [energy_dict[key][0]]
            E0.append(energy_dict[key][1].copy())#raw energy for ensemble    
            name,details,site = key.split('-') #split key into name/site
#            if 'slab' not in name: #do not include empty site energy (0)
            if site == '4_layer': #or site == '2_layer_runs':
                E0[0] -= energy_dict['slab-slab-'+site][0]
                E0[1] -= energy_dict['slab-slab-'+site][1]#subtract slab energy if adsorbed
            formula = name
            #get the composition as a list of atomic species
            try:
                composition = string2symbols(formula)
             #for each atomic species, subtract off the reference energy
                for atom in composition:
                    E0[0] = E0[0] -ref_dict[atom][0]
                    E0[1] = E0[1] -ref_dict[atom][1]
                    if site == '4_layer' and atom =='H':
                        E0[0] = E0[0] #+1.6
                        E0[1] = E0[1] #+1.6
                formation_energies[key] = E0
            except:
                if 'slab' in name:
                   formation_energies[key] = E0
            del E0
#    formation_energies['e'] = [-pot,np.zeros(2000)-pot]
    return formation_energies

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


def plot_FED(states,barrier,labels,potential):
    #make_input_file(file_name,formation_energies,frequency_dict)
    NO_rxn = MechanismPlot(states,barriers=barrier,labels=labels)
    NO_rxn.energy_mode = "absolute"
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #NO_rxn.barriers = [0]*(len(NO_rxn.energies)-1)
    NO_rxn.draw(ax=ax)
    plt.xlabel('Reaction Coordinate')
    plt.ylabel('State Free Energy')
    ax.text(0.2,0.5,'potential='+str(potential)+'V')
    fig.savefig('free_energy_diagram.jpg')
    plt.show()
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(states[4][1])
    ax.set_title("N2O2 Energy Distribution")
    fig.savefig('N2O2_hist.jpg')
    
def Surface_diagram(energy_dict,ref_offsets,surface_area,pot):
    STP_N = ref_offsets['N'][0]
    STP_N_ens = ref_offsets['N'][1]
    mu_o = STP_N - 8.617330350*10**(-5)*300*np.log(101325/101325) 
    init_mu = -278.3
    fin_mu =  -277.2
    legend_text = []
    ens_int = STP_N_ens-STP_N
    ref_offsets['N'] = [init_mu,ens_int+init_mu]
    formation_energies_ref_mv = {'N2-test-4_layer':[],
                                 'N2O-test-4_layer':[],
                                 'N2O2-test-4_layer':[],
                                 'NO-test-4_layer':[],
#                                 'H-ads_Ti-4_layer':[],
#                                 'H-ads_O-4_layer':[],
                                 'OH-test-4_layer':[],
                                 'O-test-4_layer':[],
                                 'H2O-test-4_layer':[],
#                                 'HONNOH-test-4_layer':[],
#                                 'NO2-test-4_layer':[],
                                 'slab-slab-4_layer':[]
                                 
                                }
    
    formation_energies_ens = {'N2-test-4_layer':[],
                                 'N2O-test-4_layer':[],
                                 'N2O2-test-4_layer':[],
                                 'NO-test-4_layer':[],
                                 'O-test-4_layer':[],
#                                 'H-ads_Ti-4_layer':[],
#                                 'H-ads_O-4_layer':[],
                                 'OH-test-4_layer':[],
                                 'H2O-test-4_layer':[],
#                                 'HONNOH-test-4_layer':[],
#                                 'NO2-test-4_layer':[],
                                 'slab-slab-4_layer':[]
                                 
                                }
    for i in np.arange(init_mu,fin_mu,0.01):
        formation_energies = get_formation_energies(energy_dict.copy(),ref_offsets.copy(),pot)
        for item in formation_energies_ref_mv.keys():
            formation_energies_ref_mv[item].append(formation_energies[item][0]/surface_area)
            formation_energies_ens[item].append(formation_energies[item][1].copy()/surface_area)
        ref_offsets['N'][0] = ref_offsets['N'][0] + 0.01
        ref_offsets['N'][1] = ref_offsets['N'][1] + 0.01
    
    #plotting portion
    plotting_range = np.arange(init_mu,fin_mu,0.01)-STP_N
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax2 = ax.twiny()
    ax2.set_xscale('log')
    fig.subplots_adjust(bottom=0.2)
    P = []
    for i in np.arange(init_mu,fin_mu,0.01):
        P.append(np.exp((i-mu_o)/(8.617330350*10**(-5)*300))) #Pa
    ax2.set_xlim(P[0],P[-1])
    for item in formation_energies_ref_mv.keys():
        ax2.plot(P,formation_energies_ref_mv[item])
        ax.plot(plotting_range,formation_energies_ref_mv[item])
        name,detail,site = item.split('-')
        if site =='4_layer':
            legend_text.append(name+'*')
        else:
            legend_text.append(name)
    plt.legend(legend_text)
    plt.xlabel('$\Delta \mu_N$(eV)')
    plt.ylabel('Surface Free Energy (eV/$\AA^2$)')
    ax.set_title("Surface Energy Diagram TiO$_2$")
    ax.set_xlim(plotting_range[0],plotting_range[-1])
    
    #ax2 = ax.twiny()
    ax2.set_xscale('log')
    """
    for item in formation_energies_ref_mv.keys():
        ax2.plot(P,formation_energies_ref_mv[item])
        """
    ax2.xaxis.set_ticks_position("bottom")
    ax2.xaxis.set_label_position("bottom")
    ax2.spines["bottom"].set_position(("axes", -0.15))
    plt.legend(legend_text)
    ax2.set_frame_on(True)
    ax2.patch.set_visible(False)
    for sp in ax2.spines.itervalues():
        sp.set_visible(False)
    ax2.spines["bottom"].set_visible(True)
    plt.xlabel('Pressure (atm)')
    ax.set_ylabel('Surface Free Energy (eV/$\AA^2$)')
    fig.savefig('Surface_diagram.jpg')
    
    
    #adding up the total prob for each different state
    
    
    del fig
    fig2 = plt.figure()
    ax = fig2.add_subplot(111)
    ax2 = ax.twiny()
    ax2.set_xscale('log')
    fig2.subplots_adjust(bottom=0.2)
    """for item in formation_energies_ref_mv:
        prob[item]=[]
    for item in formation_energies_ref_mv:
        for i in range(len(np.arange(init_mu,fin_mu,0.01))):
            prob[item].append(np.exp(-formation_energies_ref_mv[item][i]*surface_area/(8.617330350*10**(-5)*300)))
            for j in range(2000):
                prob[item][-1] = prob[item][-1]+np.exp(-formation_energies_ens[item][i][j]*surface_area/(8.617330350*10**(-5)*300))
    tot_P = []
    for i in range(len(np.arange(init_mu,fin_mu,0.01))):  
        tot = 0
        for item in formation_energies_ref_mv:
            tot = tot+ prob[item][i]
        tot_P.append(tot)
    for i in range(len(np.arange(init_mu,fin_mu,0.01))):   
        for item in formation_energies_ref_mv:
            prob[item][i] = prob[item][i]/tot_P[i]"""
    P_values = {}
    for key in formation_energies_ref_mv.keys():
        P_values[key]=[]
    for index in range(len(np.arange(init_mu,fin_mu,0.01))):
        P_dict = probabilities(formation_energies_ref_mv,formation_energies_ens,index)
        for key in P_dict.keys():
            P_values[key].append(P_dict[key])
    for item in formation_energies_ref_mv.keys():
        #print(plotting_range)
        #print(P_values)
        ax.plot(plotting_range,P_values[item])
        ax2.plot(P,P_values[item])
        name,detail,site = item.split('-')
        legend_text.append(name)
    plt.legend(legend_text)
    plt.xlabel('$\Delta \mu_N$(eV)')
    plt.ylabel('Probability')
    ax.set_title("Surface Probability Diagram TiO$_2$")
    ax.set_xlim(plotting_range[0],plotting_range[-1])
    ax.set_xlim(plotting_range[0],plotting_range[-1])
    #ax.set_ylim(0,1)
    ax2.set_xlim(P[0],P[-1])
    #ax2.set_ylim(0,1)
    ax2.xaxis.set_ticks_position("bottom")
    ax2.xaxis.set_label_position("bottom")
    ax2.spines["bottom"].set_position(("axes", -0.15))
    plt.legend(legend_text,loc=2)
    ax2.set_frame_on(True)
    ax2.patch.set_visible(False)
    for sp in ax2.spines.itervalues():
        sp.set_visible(False)
    ax2.spines["bottom"].set_visible(True)
    ax2.set_xlabel('Pressure (atm)')
    ax.set_xlabel('$\Delta \mu_N$(eV)')
    ax.set_ylabel('Probability')
    fig2.savefig('Prob_diagram.jpg')


P_water = vap_press('H2O',300)
#P_N2 = 0.0118168 #atm
#P_O2 = 0.0068092 #atm
energy_dict = better_pickle_reader.pickle_pull()
lem = energy_dict.copy()
#energy_dict = better_pickle_reader.add_ZPE(energy_dict)
#print(energy_dict['NO2-molecules-DFT_references'][0])
#print(energy_dict['H2O-molecules-DFT_references'][0])
energy_dict = better_pickle_reader.add_vib_energy(energy_dict,300,101325,
                                                  spec_case=[
                                                  'vib-H2O-molecules-DFT_references',
#                                                  'vib-N2-molecules-DFT_references',
#                                                  'vib-O2-molecules-DFT_references'
                                                  ],
                                                 spec_press=[P_water*101325,
#                                                              P_N2*101324,
#                                                              P_O2*101325
                                                              ]
                                                  )
                                                  
print(energy_dict['N2-test-4_layer'][0])
#print (2*energy_dict['O-test-4_layer'][0]-energy_dict['O2-molecules-DFT_references'][0]-2*energy_dict['slab-slab-4_layer'][0])
#print (energy_dict['O2-test-4_layer'][0]-energy_dict['O2-molecules-DFT_references'][0]-energy_dict['slab-slab-4_layer'][0])
#print (energy_dict['N2O2-test-4_layer'][0]-2*energy_dict['NO-molecules-DFT_references'][0]-energy_dict['slab-slab-4_layer'][0])
#print (energy_dict['NO-test-4_layer'][0]-energy_dict['NO-molecules-DFT_references'][0]-energy_dict['slab-slab-4_layer'][0])
#print (energy_dict['N2O-test-4_layer'][0]-energy_dict['N2O-molecules-DFT_references'][0]-energy_dict['slab-slab-4_layer'][0])
#print(2*energy_dict['NO-molecules-DFT_references'][0]-energy_dict['O2-molecules-DFT_references'][0]-energy_dict['N2-molecules-DFT_references'][0])

"""
print(energy_dict['NO-molecules-DFT_references'][0]+energy_dict['H2-molecules-DFT_references'][0]-energy_dict['H2O-molecules-DFT_references'][0]-0.5*energy_dict['N2-molecules-DFT_references'][0])
print(energy_dict['H2O-molecules-DFT_references'][0]-0.5*energy_dict['O2-molecules-DFT_references'][0]-energy_dict['H2-molecules-DFT_references'][0])
print(energy_dict['NO-molecules-DFT_references'][0]-0.5*energy_dict['O2-molecules-DFT_references'][0]-0.5*energy_dict['N2-molecules-DFT_references'][0])
l = energy_dict['H2O-molecules-DFT_references'][1]-0.5*energy_dict['O2-molecules-DFT_references'][1]-energy_dict['H2-molecules-DFT_references'][1]
h = energy_dict['NO-molecules-DFT_references'][1]-0.5*energy_dict['O2-molecules-DFT_references'][1]-0.5*energy_dict['N2-molecules-DFT_references'][1]
print(np.std(-l+h)/2)
print(np.average(-l+h)/2)"""
atomic_references = {'N':'N2','O':'O2','H':'H2O','C':'CO2'}
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
    if site == 'DFT_references' and 'molecules_2' not in key:
        compound_en_dict[name] = [energy_dict[key][0]]
        compound_en_dict[name].append(energy_dict[key][1].copy())
"""
initial_state = ['slab-slab-4_layer','N2-molecules-DFT_references','O2-molecules-DFT_references']
states = [['N2-test-4_layer','O2-molecules-DFT_references'],
          ['N2O-test-4_layer','0.5*O2-molecules-DFT_references'],
          ['N2O2-test-4_layer'],
          ['2*NO-test-4_layer'],
          ['2*NO-molecules-DFT_references']
          ] 
"""
initial_state = ['slab-slab-4_layer','2*H2O-molecules-DFT_references']
states = [['N2-test-4_layer','2*H2O-molecules-DFT_references'],
          ['N2-test-4_layer','OH-test-4_layer'],
          ['N2H-test-4_layer','H2O-molecules-DFT_references'],#'H2-molecules-DFT_references'#,'2*e'],
          ['HNNH-test-4_layer'],#'2*H2-molecules-DFT_references'],#,'4*e'],
          ['H2NNH2-test-4_layer'],#'2*H2-molecules-DFT_references'],#,'4*e'],
          #['H2NNH2-test-4_layer'],#'2*H2-molecules-DFT_references']#,'4*e']
          ['NH3-test-4_layer','NO-test-4_layer'],#'2*H2-molecules-DFT_references']#,'4*e']
          #['NH3-test-4_layer','NO-test-4_layer'],#'2*H2-molecules-DFT_references']#,'4*e']
          ]
          
pot = 0
ref_test, ref_offsets = reference_maker.convert_formation_energies(compound_en_dict.copy(),atomic_references,composition_dict,pot)
#print ref_offsets['H'][0]
#print ref_offsets['O'][0]
#print formation_energies_ref_mv.keys()

formation_energies_ref_mv = {'N2-test-4_layer':[],
                                 'N2O-test-4_layer':[],
                                 'N2O2-test-4_layer':[],
                                 'NO-test-4_layer':[],
#                                 'OHNNHO-ads_Ti-4_layer':[],
                                 'O-test-4_layer':[],
#                                 'H-ads_O-4_layer':[],
                                 'OH-test-4_layer':[],
                                 'H2O-test-4_layer':[]
                                 }
                                 
formation_energies_ref_l =     {'N2-test-4_layer':[],
                                 'N2O-test-4_layer':[],
                                 'N2O2-test-4_layer':[],
                                 'NO-test-4_layer':[],
#                                 'OHNNHO-ads_O-4_layer':[],
#                                 'H-ads_O-4_layer':[],
                                 'OH-test-4_layer':[],
                                 'H2O-test-4_layer':[],
                                 'O-test-4_layer':[],
                                 'HNNOH-test-4_layer':[],
#                                 'NO-2_NO-4_layer':[],
                                 'slab-slab-4_layer':[]
                                    }
formation_energies = get_formation_energies(energy_dict.copy(),ref_offsets.copy(),pot)
#print('species\t relative energy\t stdev')
for item in formation_energies.keys():
    if abs(formation_energies[item][0])<50:
        print('\''+item+'\'\t'+str(formation_energies[item][0])+'\t'+str(np.std(formation_energies[item][1])))
        #print(formation_energies[item][0])
        n=1

#print(np.mean(formation_energies['HONNOH-test-4_layer'][1]))
#print formation_energies['N2-test-4_layer']
for item in formation_energies_ref_mv:
    formation_energies_ref_mv[item] = formation_energies[item]

#print(formation_energies['H2O-molecules-DFT_references'][0]-formation_energies['H2-molecules-DFT_references'][0]-0.5*formation_energies['O2-molecules-DFT_references'][0])
#Gibbs_isotherm(formation_energies_ref_mv.copy(),'N2-test-4_layer')
#print formation_energies['H2-molecules-DFT_references']
states = rxn_energy(states,initial_state,formation_energies.copy())
#print(states[0][0]-states[-1][0])

initial_state2 = ['slab-slab-4_layer','N2-molecules-DFT_references','3*H2-molecules-DFT_references']
states2 =[
          ['N2-test-4_layer','3*H2-molecules-DFT_references'],
          ['N2H-test-4_layer','2.5*H2-molecules-DFT_references'],
          ['HNNH-test-4_layer','2*H2-molecules-DFT_references'],
          ['H2NNH-test-4_layer','1.5*H2-molecules-DFT_references'],
          ['H2NNH2-test-4_layer','NH2-test-4_layer','0.5*H2-molecules-DFT_references'],
          ['2*NH3-test-4_layer'],
          ['slab-slab-4_layer','NH3-molecules-DFT_references'] 
          ]
states2 = rxn_energy(states2,initial_state2,formation_energies)


plot_FED(states,[-15,-15,-15,-15,-15],['N2+O2', 'N2*+O2','N2O*+0.5O2','N2O2*','2NO*','2NO+*'],pot)
#plot_FED(states2,[-15*len(states2)],['N2+3H2',''])
#ref_test2, ref_offsets2 = reference_maker.convert_formation_energies(compound_en_dict.copy(),atomic_references,composition_dict,potential=0)
surface_area = 6.5797272*5.96159647
Surface_diagram(energy_dict.copy(),ref_offsets.copy(),surface_area,pot)

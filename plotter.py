# -*- coding: utf-8 -*-
"""
Created on Fri May  5 16:03:35 2017

@author: benjamin
"""

import sys
from ase.atoms import string2symbols
import reader
import reference_maker
from matplotlib import cm
import matplotlib.ticker as ticker
from Pourbaix import oxidation_state
from catmap.analyze import MechanismPlot
#import catmap.data.parameter_data#.experimental_gas_formation_energies as egfe
import numpy as np
from mapping import mapping_dict
import pylab as plt
from ThemoChem import vap_press
from ase.units import kB
from Beef_data import beef as bf
from reference_maker import get_formation_energies
cmap = cm.Dark2
plt.locator_params(axis='y', nticks=1)
plt.locator_params(axis='x', nticks=1)
from pylab import MaxNLocator


#plt = catmap.plt

#P_N2 = 0.0118168 #atm
#P_O2 = 0.0068092 #atm
plt.rc('font', size=10)
plt.rc('axes', titlesize=12)
plt.rc('xtick', labelsize=12)    # fontsize of the tick labels
plt.rc('ytick', labelsize=12)    # fontsize of the tick labels
plt.rc('legend', fontsize=10)
#font = {'family' : 'normal',
#        'size'   : 6}
font = {}
figsize_dict = {'figsize':(4,3)}
#plt.rc('font', **font)
cs = ['r','g','b','k','m','c',cmap(0.5),cmap(0.6),cmap(0.7)]
lin = 'k'
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


def plot_FED(states,barrier,labels,potential,name,atomic_references):
    #make_input_file(file_name,formation_energies,frequency_dict)
    NO_rxn = MechanismPlot(states,barriers=barrier,labels=labels)
    NO_rxn.energy_mode = "absolute"
    
    fig = plt.figure(**figsize_dict)
    ax = fig.add_subplot(111)
    #NO_rxn.barriers = [0]*(len(NO_rxn.energies)-1)
    NO_rxn.draw(ax=ax)
    plt.xlabel('Reaction Coordinate')
    plt.ylabel('State Free Energy (eV)',labelpad=0)
    textstr='potential='+str(potential)+'V\n'#+'References:\n'
#    for item in atomic_references.keys():
#        if item !='C':
#            textstr+=item+':'+atomic_references[item]+'\n'
    #ax.text(0.2,-0.2,'potential='+str(potential)+'V')
    props = dict(facecolor='white', alpha=0.5)
    ax.text(0.05, 0.08, textstr, transform=ax.transAxes, fontsize=14,
        verticalalignment='top')
    ya = ax.get_yaxis()
    ya.set_major_locator(MaxNLocator(integer=True))
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.xaxis.set_minor_locator(ticker.NullLocator())


    fig.savefig(name+'.pdf')
    plt.show()
    
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    ax.hist(states[4][1])
#    ax.set_title("N2O2 Energy Distribution")
#    fig.savefig('N2O2_hist.jpg')
    #return fig.axes[0]
    
    
    
def plot_FED_return(states,barrier,labels,potential,atomic_references,figure_obj,share=None):
    #make_input_file(file_name,formation_energies,frequency_dict)
    NO_rxn = MechanismPlot(states,barriers=barrier,labels=labels)
    NO_rxn.energy_mode = "absolute"
    
    #ax = figure_obj.add_subplot(inp,sharex=share)
    #NO_rxn.barriers = [0]*(len(NO_rxn.energies)-1)
    NO_rxn.draw(ax=figure_obj)
    plt.xlabel('Reaction Coordinate')
    plt.ylabel('State Free Energy (eV)')

    textstr='potential='+str(potential)+'V\n'+'References:\n'
    for item in atomic_references.keys():
        if item !='C':
            textstr+=item+':'+atomic_references[item]+'\n'
    #ax.text(0.2,-0.2,'potential='+str(potential)+'V')
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    figure_obj.text(0.05, 0.05, textstr, transform=figure_obj.transAxes, fontsize=14,
        verticalalignment='bottom', bbox=props)
    figure_obj.xaxis.set_major_locator(ticker.NullLocator())
    figure_obj.xaxis.set_minor_locator(ticker.NullLocator())
    return figure_obj

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

def Gibbs_isotherm(energy_dict,ens_dict,index):
    k = 8.617e-5
    T = 300
    a = 6.5797272*5.96159647
    for key in energy_dict.keys():
        new_entry = ens_dict[key][index]
        #new_entry = np.append(new_entry,energy_dict[key][index])
        energy_dict[key] = new_entry
    P_dict = {}
    for sp,ens in energy_dict.items():
        P_dict[sp] = []
        e = np.mean(ens)
        x = np.sum(np.exp(-e*a/(k*T)))
        P_dict[sp] += [x]

    full_dict = {}
    for sp in P_dict:
        full_dict[sp] = 0
        for i in range(len(P_dict[sp])):
            full_dict[sp] += P_dict[sp][i]/sum([P_dict[sp_i][i] for sp_i in P_dict])
        full_dict[sp] /= len(P_dict[sp]) #normalize
    return full_dict

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

def get_formation_energies_photochemistry(energy_dict1,ref_dict,band_details):
    """
    band details should be a list [band_gap,upper_edge]
    """
    ue = band_details[1]
    le = (band_details[1]-band_details[0])
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
                oxy = oxidation_state(formula)
                #print(key)
                EU = E0[0] - oxy*ue
                EL = E0[0] - oxy*le
                if EU == EL:
                    pot = 0
                elif EU<EL:
                    pot = ue
                else:
                    pot = le
                E0[0] = E0[0] - oxy*pot
                E0[1] = E0[1] - oxy*pot
             #for each atomic species, subtract off the reference energy
                for atom in composition:
                    E0[0] = E0[0] -ref_dict[atom][0]
                    E0[1] = E0[1] -ref_dict[atom][1]
                    if site == '4_layer' and atom =='H':
                        E0[0] = E0[0] 
                        E0[1] = E0[1] 
                formation_energies[key] = E0 #+oxy*pot
            except:
                if 'slab' in name:
                   formation_energies[key] = E0
            del E0
    #formation_energies['e'] = [-pot,np.zeros(2000)-pot]
    return formation_energies
    
def move_free_energy(energy_dict,ref_offsets,surface_area,atom_to_move,relative_molecule):
    pot=0
    step = 0.001
    ll= 0.55 #lower limit
    ul = 0.2 #upper limit
    mv = atom_to_move
    mult,atm = mv.split('*')
    mult = float(mult)
    mol = string2symbols(relative_molecule)
    STP = 0
    STP_ens = 0
    for atom in mol:
        STP += ref_offsets[atom][0]
        STP_ens += ref_offsets[atom][1]
#    STP = ref_offsets[atm][0]
#    STP_ens = mult*ref_offsets[atm][1]
    mu_o = STP - 8.617330350*10**(-5)*300*np.log(101325/101325)
    init_mu_cmp = STP-ll
    fin_mu_cmp =  STP+ul
    init_mu = STP - (mult*ll)
    fin_mu =  STP + (mult*ul)
    ens_int = STP_ens-STP
    ro = ref_offsets.copy()
    ref_offsets[atm] = [ro[atm][0]-mult*ll,ro[atm][1]-mult*ll]
    formation_energies_ref_mv = {
                                 'N2-test-4_layer':[],
#                                 'N2O-test-4_layer':[],
#                                 'N2O2-test-4_layer':[],
#                                 'NO-test-4_layer':[],
#                                 'H-ads_Ti-4_layer':[],
#                                 'H-ads_O-4_layer':[],
                                 'OH-test-4_layer':[],
                                 'H2O-test-4_layer':[],
                                 'slab-slab-4_layer':[],
                                 'H4O2-2_H2O-4_layer':[],
                                 'N4-N2_1ML-4_layer':[],
                                 'O4-full_O2_coverage-4_layer':[],
                                 'O2-test-4_layer':[],  
                                 
                                }
#    
    formation_energies_ens = {
                                 'N2-test-4_layer':[],
#                                 'N2O-test-4_layer':[],
#                                 'N2O2-test-4_layer':[],
#                                 'NO-test-4_layer':[],
#                                 'H-ads_Ti-4_layer':[],
#                                 'H-ads_O-4_layer':[],
                                 'OH-test-4_layer':[],
                                 'H2O-test-4_layer':[],
                                 'slab-slab-4_layer':[],
                                 'H4O2-2_H2O-4_layer':[],
                                 'N4-N2_1ML-4_layer':[],
                                 'O4-full_O2_coverage-4_layer':[],
                                 'O2-test-4_layer':[],             
                                }
#    for item in energy_dict:
#        formation_energies_ref_mv[item]=[]
#        formation_energies_ens[item]=[]
    for i in np.arange(init_mu,fin_mu,step*mult):
        formation_energies = get_formation_energies(energy_dict.copy(),ref_offsets.copy(),pot)
        for item in formation_energies_ref_mv.keys():
            formation_energies_ref_mv[item].append(formation_energies[item][0]/surface_area)
            formation_energies_ens[item].append(formation_energies[item][1].copy()/surface_area)
        ref_offsets[atm][0] = ref_offsets[atm][0] + (step)*mult
        ref_offsets[atm][1] = ref_offsets[atm][1] + (step)*mult
    P = []
    for i in np.arange(init_mu_cmp,fin_mu_cmp,step):
        P.append(np.exp((i-mu_o)/(8.617330350*10**(-5)*300))) #Pa
    return formation_energies_ref_mv,formation_energies_ens,STP,np.arange(init_mu_cmp,fin_mu_cmp,step),P

def plot_surface_energy_diagram(formation_energies_ref_mv,P,cmp_mu_range,STP,cmpd):    
    #plotting portion
    legend_text =[]
    plotting_range = cmp_mu_range
    fig = plt.figure(**figsize_dict)
    ax = fig.add_subplot(111)
    ax2 = ax.twiny()
    ax2.set_xscale('log')
    fig.subplots_adjust(bottom=0.2)
    ax2.set_xlim(P[0],P[-1])
    d = len(formation_energies_ref_mv.keys())
    N=0
    for item in sorted(formation_energies_ref_mv.keys()):
        ax2.plot(P,formation_energies_ref_mv[item],color=cs[N])
        ax.plot(plotting_range,formation_energies_ref_mv[item],color=cs[N])
        name,detail,site = item.split('-')
        N=N+1
#        if site =='4_layer' and name!= 'slab':
#            legend_text.append(name+'*')
#        else:
#            legend_text.append(name)
        if site =='4_layer':
            legend_text.append(mapping_dict[item])
    plt.legend(legend_text)
    ax.set_xlabel('$\Delta \mu_{'+cmpd+'}$(eV)',labelpad=0)
    plt.ylabel('Surface Free Energy (eV/$\AA^2$)')
    ax.set_title("Surface Energy Diagram TiO$_2$")
    ax.set_xlim(plotting_range[0]-STP,plotting_range[-1]-STP)
    
    if cmpd == 'H_2O':
        ax2.plot([0.0357/1.01325,0.0357/1.01325],ax2.get_ylim(),color=lin,ls='--')
        mu = kB*300*np.log(0.0357/1.01325)
        ax.plot([mu,mu],ax.get_ylim(),color=lin,ls='--')
        ax.text(mu, 0.85, '$H_2O_{(l)}$')
    elif cmpd =='N_2':
        ax2.plot([0.0118168,0.0118168],ax2.get_ylim(),color=lin,ls='--')
        mu = kB*300*np.log(0.0118168)
        ax.plot([mu,mu],ax.get_ylim(),color=lin,ls='--')
        ax.text(mu, -0.02, 'Henry\'s Law')
    #ax2 = ax.twiny()
    ax2.set_xscale('log')
    """
    for item in formation_energies_ref_mv.keys():
        ax2.plot(P,formation_energies_ref_mv[item])
        """
    ax2.xaxis.set_ticks_position("bottom")
    ax2.xaxis.set_label_position("bottom")
    ax2.spines["bottom"].set_position(("axes", -0.15))
    ## tick locations
    ax.xaxis.set_ticks([-0.4,-0.2,0,0.2])
    ax2.xaxis.set_ticks([10**-9,10**-5,10**-1,10**3])
    ax2.yaxis.set_ticks([0.12,0.08,0.04,0,-0.04])
    ax.yaxis.set_ticks([0.12,0.08,0.04,0,-0.04])
    
    plt.legend(legend_text)
    ax2.set_frame_on(True)
    ax2.patch.set_visible(False)
    for sp in ax2.spines.itervalues():
        sp.set_visible(False)
    ax2.spines["bottom"].set_visible(True)
    plt.xlabel('Pressure (atm)',labelpad=0)
    ax.set_ylabel('Surface Free Energy (eV/$\AA^2$)',labelpad=0)
    return fig
    
def Probability(formation_energies_ref_mv, formation_energies_ens, cmp_mu_range,P,surface_area,STP,cmpd):
    legend_text = []
    fig2 = plt.figure(**figsize_dict)
    ax = fig2.add_subplot(111)
    ax2 = ax.twiny()
    ax2.set_xscale('log')
    fig2.subplots_adjust(bottom=0.2)
    P_values = {}
    d = len(formation_energies_ref_mv.keys())
    N=-1
    for key in formation_energies_ref_mv.keys():
        P_values[key]=[]
    for index in range(len(cmp_mu_range)):
        P_dict = probabilities(formation_energies_ref_mv,formation_energies_ens,index)
        for key in P_dict.keys():
            P_values[key].append(P_dict[key])
    for item in sorted(formation_energies_ref_mv.keys()):
        #print(plotting_range)
        #print(P_values)
        N = N+1
        ax.plot(cmp_mu_range,P_values[item],color=cs[N])
        ax2.plot(P,P_values[item],color=cs[N])
        name,detail,site = item.split('-')
#        if site =='4_layer' and name!= 'slab':
#            legend_text.append(name+'*')
#        else:
#            legend_text.append(name)
        if site =='4_layer':
            legend_text.append(mapping_dict[item])
    if cmpd == 'H_2O':
        ax2.plot([0.0357/1.01325,0.0357/1.01325],ax2.get_ylim(),color=lin,ls='--')
        mu = kB*300*np.log(0.0357/1.01325)
        ax.plot([mu,mu],ax.get_ylim(),color=lin,ls='--')
        ax.text(mu, 0.85, '$H_2O_{(l)}$')
    elif cmpd =='N_2':
        ax2.plot([0.0118168,0.0118168],ax2.get_ylim(),color=lin,ls='--')
        mu = kB*300*np.log(0.0118168)
        ax.plot([mu,mu],ax.get_ylim(),color=lin,ls='--')
        ax.text(mu, 0.85, 'Henry\'s Law')   
    plt.legend(legend_text)
    plt.xlabel('$\Delta \mu_{'+cmpd+'}$(eV)')
    plt.ylabel('Probability')
    ax.set_title("Surface Probability Diagram TiO$_2$")
    ax.set_xlim(cmp_mu_range[0]-STP,cmp_mu_range[-1]-STP)
    ax.set_xlim(cmp_mu_range[0]-STP,cmp_mu_range[-1]-STP)
    ## tick locations
    ax.xaxis.set_ticks([-0.4,-0.2,0,0.2])
    ax2.xaxis.set_ticks([10**-9,10**-5,10**-1,10**3])
    ax2.yaxis.set_ticks([1,0.75,0.5,0.25,0])
    ax.yaxis.set_ticks([1,0.75,0.5,0.25,0])
    
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
    ax2.set_xlabel('Pressure (atm)',labelpad=0)
    ax.set_xlabel('$\Delta \mu_{'+cmpd+'}$(eV)',labelpad=0)
    ax.set_ylabel('Probability')
    #fig2.savefig('Prob_diagram.jpg')
    return fig2

def Coverage(formation_energies_ref_mv, formation_energies_ens, cmp_mu_range,P,surface_area,STP,cmpd):
    legend_text = []
    fig2 = plt.figure(**figsize_dict)
    ax = fig2.add_subplot(111)
    ax2 = ax.twiny()
    ax2.set_xscale('log')
    fig2.subplots_adjust(bottom=0.2)
    P_values = {}
    for key in sorted(formation_energies_ref_mv.keys()):
        P_values[key]=[]
    for index in range(len(cmp_mu_range)):
        P_dict = Gibbs_isotherm(formation_energies_ref_mv,formation_energies_ens,index)
        for key in P_dict.keys():
            P_values[key].append(P_dict[key])
    d = len(formation_energies_ref_mv.keys())
    N=0
    for item in sorted(formation_energies_ref_mv.keys()):
        #print(plotting_range)
        #print(P_values)
        #print(c)
        ax.plot(cmp_mu_range,P_values[item],color=cs[N])
        ax2.plot(P,P_values[item],color=cs[N])
        name,detail,site = item.split('-')
        N = N+1
#        if site =='4_layer' and name!= 'slab':
#            legend_text.append(name+'*')
#        else:
#            legend_text.append(name)
        if site =='4_layer':
            legend_text.append(mapping_dict[item])
    if cmpd == 'H_2O':
        ax2.plot([0.0357/1.01325,0.0357/1.01325],ax2.get_ylim(),color=lin,ls='--')
        mu = kB*300*np.log(0.0357/1.01325)
        ax.plot([mu,mu],ax.get_ylim(),color=lin,ls='--')
        ax.text(mu, 0.85, '$H_2O_{(l)}$')
    elif cmpd =='N_2':
        ax2.plot([0.0118168,0.0118168],ax2.get_ylim(),color=lin,ls='--')
        mu = kB*300*np.log(0.0118168)
        ax.plot([mu,mu],ax.get_ylim(),color=lin,ls='--')
        ax.text(mu, 0.85, 'Henry\'s Law')
    plt.legend(legend_text)
    plt.xlabel('$\Delta \mu_{'+cmpd+'}$(eV)',labelpad=0)
    plt.ylabel('$\Theta$ Coverage',labelpad=0)
    ax.set_title("Surface Coverage Diagram TiO$_2$")
    ax.set_xlim(cmp_mu_range[0]-STP,cmp_mu_range[-1]-STP)
    ax.set_xlim(cmp_mu_range[0]-STP,cmp_mu_range[-1]-STP)
     #tick locations  
    ax.xaxis.set_ticks([-0.4,-0.2,0,0.2])
    ax2.xaxis.set_ticks([10**-9,10**-5,10**-1,10**3])
    ax2.yaxis.set_ticks([1,0.75,0.5,0.25,0])
    ax.yaxis.set_ticks([1,0.75,0.5,0.25,0])
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
    ax.set_xlabel('$\Delta \mu_{'+cmpd+'}$(eV)',labelpad=0)
    ax.set_ylabel('$\Theta$ Coverage',labelpad=0)
    #fig2.savefig('coverage_diagram.jpg')
    return fig2

def strip_high_energy(eng,ens,cutoff): 
    #strips out species that are higher than the cutoff over the entire range
    for item in eng.keys():
        if 'DFT' in item:
            del eng[item]
            del ens[item]
    random_key = eng.keys()[0]
    lowest_e = [0]*len(eng[random_key])
    TF_dict = {}
    for item in eng:
        TF_dict[item] = [0]*len(eng[random_key])
    for i in range(len(eng[random_key])):
        lowest_e[i]=eng[random_key][i]
        for item in eng:
            if eng[item][i]<lowest_e[i]:
                lowest_e[i]=eng[item][i]
#    fig2 = plt.figure()
#    ax = fig2.add_subplot(111)
#    ax.plot(range(len(eng[random_key])),lowest_e)
#    fig2.savefig('lol.jpg')
    
    for i in range(len(eng[random_key])):
        for item in eng.keys():
            if abs(eng[item][i]-lowest_e[i])>cutoff:
                TF_dict[item][i]=0
            else:
                TF_dict[item][i]=1
    
    for item in eng.keys():
        if 1 not in TF_dict[item]:
            del eng[item]
            del ens[item]
    return eng,ens

def ref_to_OH_hack(ref_offsets,energy_dict):
    ed = energy_dict.copy()
    ref_offsets['H'][0]=0.5*(5.63-ed['H2-molecules-DFT_references'][0]+2*ed['H2O-molecules-DFT_references'][0])-0.5*ed['O2-molecules-DFT_references'][0]
    ref_offsets['H'][1]=0.5*(5.63-ed['H2-molecules-DFT_references'][1]+2*ed['H2O-molecules-DFT_references'][1])-0.5*ed['O2-molecules-DFT_references'][1]
    return ref_offsets
    
def eng_add_OH_hack(energy_dict):
    ed = energy_dict.copy()
    eng = 0.5*(5.4521-ed['H2-molecules-DFT_references'][0]+2*ed['H2O-molecules-DFT_references'][0])
    ens = 0.5*(5.4521-ed['H2-molecules-DFT_references'][1]+2*ed['H2O-molecules-DFT_references'][1])
    energy_dict['OH-H2_H2O_hacked-DFT_references'] = bf([eng,ens])
#    energy_dict['OH-H2_H2O_hacked-DFT_references'] = []
#    energy_dict['OH-H2_H2O_hacked-DFT_references'].append(0.5*(5.4521-ed['H2-molecules-DFT_references'][0]+2*ed['H2O-molecules-DFT_references'][0]))
#    energy_dict['OH-H2_H2O_hacked-DFT_references'].append(0.5*(5.4521-ed['H2-molecules-DFT_references'][1]+2*ed['H2O-molecules-DFT_references'][1]))
    return energy_dict
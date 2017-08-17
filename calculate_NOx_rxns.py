import numpy as np
from glob import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

### DEFINE FUNCTIONS

def get_BEEF_ensemble(log):
    f = open(log)
    txt = f.read()
    f.close()
    _,E_total = txt.rsplit('total energy              =',1)
    E_total,_ = E_total.split('Ry',1)
    E_total = float(E_total.strip())
    E_total *= 13.605698

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

### DEFINE DATA

nist = {
    'CH3OCH3': -3.831687734, 
    'isobutene': -2.7612012400000001, 
    'CH3CH2OH': -4.3257274452000001, 
    'H2CCO': -1.2932338992999999, 
    'CH3COOH': -5.936414407, 
    'butadiene': -0.94117440899999982, 
    'CH3NO2': -1.9433241520000002, 
    'CH2NHCH2': -0.37361471699999993, 
    'C3H8': -3.5838955599999998, 
    'CH3OH': -3.3154954364, 
    'CH2OCH2': -0.70501341749999991, 
    'CCH': 5.4911384108000005, 
    'CH3CH2NH2': -2.7060065129999997, 
    'N2O': 0.59071737410000003, 
    'C3H6_Cs': -1.7302961269999999, 
    'O3': 1.3200076574999999, 
    'O2': -0.097959758399999999, 
    'NO2': 0.15269533860000001, 
    'CH2_s1A1d': 3.6179179609999998, 
    'NH': 3.5110470303999999, 
    'NCCN': 2.7613439944999998, 
    'CH3O': -0.57019516400000003, 
    'H2O2': -2.0341974058000001, 
    'cyclobutene': -0.49233188199999978, 
    'NO': 0.8203358394000001, 
    'OCHCHO': -3.1146914399999996, 
    'cyclobutane': -2.3477410291999998, 
    'CH3CN': -0.45205272699999999, 
    'C2H2': 1.6608275572000002, 
    'CH3': 0.7658170968000001, 
    'CH4': -1.8659396301, 
    'CH3COCH3': -4.2818641334, 
    'CH': 5.9634832613999995, 
    'CO': -1.3140678102999999, 
    'CN': 4.4009493388000003, 
    'isobutane': -4.5367431759999999, 
    'C4H4O': -2.0735235359999997, 
    'C3H4_C3v': 0.52048922799999997, 
    'NH3': -1.2981708845000002, 
    'NH2': 1.4628439994, 
    'CH3CHO': -3.1157422380000002, 
    'C2H6': -2.6688338305999997, 
    'N2': -0.14621433119999999, 
    'C2H4': -0.70501341749999991, 
    'HCN': 0.94892666259999991, 
    'HCO': 0.097736827300000051, 
    'CH3CONH2': -4.2117845308000001, 
    'H3CNH2': -1.7742767699999999, 
    'CO2': -4.3855609798000001, 
    'OH': 0.1529043621, 
    'H2': -0.27283919039999999, 
    'methylenecyclopropane': 0.0057401610000002989, 
    'C6H6': -1.6117550120000002, 
    'N2H4': -0.25599903979999983, 
    'C4H4NH': -0.88607863680000021, 
    'H2O': -3.0344480389999999, 
    'H2CCHCN': 0.60457370399999988, 
    'HCOOH': -4.7369155476999998, 
    'H2CO': -1.7865010562000001, 
    'graphite':-0.1596,
    }

basepath = '/nv/hp22/amedford6/medford-shared/users/bcomer3/espresso_rutile/espresso_Rutile/molecules'
allpaths = glob(os.path.join(basepath, '*', 'esp.log', 'log'))
ensembles = {}

for p in allpaths:
    name = p.split(os.sep)[-3]
    ensembles[name] = get_BEEF_ensemble(p)[-1]

### ANALYZE DATA

rxns = {}
rxns['NO'] = [['NO'],['0.5*N2','0.5*O2']]
rxns['NO2'] = [['NO2'],['0.5*N2','O2']]
rxns['N2O'] = [['N2O'],['0.5*O2','N2']]
rxns['H2O'] = [['H2O'],['0.5*O2','H2']]
rxns['CO2'] = [['CO2'],['0.5*O2','CO']]

total_error = 0
for rxn,p_r in rxns.items():
    print rxn
    p,r = p_r
    err = get_error(p,r,ensembles,nist)
    print err.mean(), err.std(), err.min()
    total_error += err

print 'Total:'
print total_error.mean(), total_error.std(), total_error.min()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(total_error)
ax.set_xlabel('Total Error [eV]')
ax.set_ylabel('Counts')
fig.savefig('total_error_histogram.pdf')

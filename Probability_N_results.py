import numpy as np
import pylab as plt
from copy import deepcopy as copy

sigma_H2O = 0.25
sigma_N2 = 0.30
sigma_N2O = 0.1
N = 2000
surface_area = 6.5797272*5.96159647
energy_dict = {  #first coulmn is relative energy, second is stdev

'H2O':[-0.108642605097,0.210366264304],
'N2O':[0.287919113635,0.197975239689],
'NO':[0.423551909117,0.202128849461],
'N2':[0.132870077634,0.185160997229],
#'N2O2':[0.806231316239,0.437702049029],
'H':[3.63761346653,0.124960371558],
'HONNOH':[1.47052737239,0.427984423178],
'slab':[0.0,0.0],
'OH':[0.459312245773,0.265970368576],
}

for item in energy_dict:
    energy_dict[item] = (energy_dict[item][0] + np.random.normal(0, energy_dict[item][1], N))/surface_area

def probabilities(energy_dict1):
    k = 8.617e-5
    T = 300
    P_dict = {}
    energy_dict = energy_dict1.copy()
    for sp,ens in energy_dict.items():
        P_dict[sp] = []
        for i, e in enumerate(ens):
            x = np.sum(np.exp(-e/(k*T)))
            P_dict[sp] += [x]
    full_dict = {}
    for sp in P_dict:
        full_dict[sp] = 0
        for i in range(len(P_dict[sp])):
            full_dict[sp] += P_dict[sp][i]/sum([P_dict[sp_i][i] for sp_i in P_dict])
        full_dict[sp] /= len(P_dict[sp]) #normalize
    return full_dict


def shift_energies(energy_dict, mu_N):
    new_dict = energy_dict.copy()
#    new_dict['N2O'] = new_dict['N2O']-2*mu_N
#    new_dict['N2O2'] = new_dict['N2O2']-2*mu_N
    new_dict['N2'] = new_dict['N2']-2*mu_N
    new_dict['NO'] = new_dict['NO']-mu_N
    new_dict['N2O'] = new_dict['N2O']-2*mu_N
    new_dict['HONNOH'] = new_dict['HONNOH']-2*mu_N
    return new_dict

mu_range = np.linspace(-1.5,2.2,100)
P_values = {}
for item in energy_dict:
    P_values[item]= []
for mu in mu_range.copy():
    new_energies = shift_energies(energy_dict.copy(), mu)
    P_dict = probabilities(new_energies)
    for key in P_dict.keys():
        P_values[key].append(P_dict[key])

fig = plt.figure()
ax = fig.add_subplot(111)

for key in P_values:
    ax.plot(mu_range, P_values[key],label=key)

plt.legend()
plt.show()
import numpy as np
import pylab as plt

sigma_H2O = 0.25
sigma_N2 = 0.30
sigma_N2O = 1
N = 2000
E_H2O = -3.5 +  np.random.normal(0, sigma_H2O, N) #simulated ensemble H2O
E_N2 = -3.45 +  np.random.normal(0, sigma_N2, N) #simulated ensemble N2
E_N2O = -3.40 +  np.random.normal(0, sigma_N2O, N) #simulated ensemble N2

energy_dict = {'N2':E_N2,'H2O':E_H2O,'N2O':E_N2O}

#for key, val in energy_dict.items():
    #print key, val.mean()

def probabilities(energy_dict):
    k = 8.617e-5
    T = 300
#    for key in energy_dict.keys():
#        new_entry = energy_dict[key][1]
#        new_entry = np.append(new_entry,energy_dict[key][0])
#        energy_dict[key] = new_entry
    P_dict = {}
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
    new_dict['N2'] = new_dict['N2']-2*mu_N
    new_dict['N2O'] = new_dict['N2O']-2*mu_N
    return new_dict

mu_range = np.linspace(-0.5,1,100)
P_values = {'N2':[], 'H2O':[],'N2O':[]}
for mu in mu_range:
    new_energies = shift_energies(energy_dict, mu)
    P_dict = probabilities(new_energies)
    for key in P_dict.keys():
        P_values[key].append(P_dict[key])

fig = plt.figure()
ax = fig.add_subplot(111)

for key in P_values:
    ax.plot(mu_range, P_values[key],label=key)
#plt.show()
plt.legend()
plt.savefig('surf.jpg')
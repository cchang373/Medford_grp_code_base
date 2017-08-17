# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 16:26:16 2017

@author: benjamin
"""
import numpy as np
#np.set_printoptions(threshold=np.nan)
log = 'log'
f = open(log)
txt = f.read()
f.close()
_,E_total = txt.rsplit('total energy              =',1)
E_total,_ = E_total.split('Ry',1)
E_total = float(E_total.strip())
E_total *= 13.605698
_, ens = txt.rsplit('BEEFens 2000 ensemble energies',1)
ens,_ = ens.split('BEEF-vdW xc energy contributions',1)
g = open('ensemble.txt', 'wr')
ens.strip()
g = open('ensemble.txt', 'wr')
g.write(ens)
g.close()
ens_ryd = []
for Ei in ens.split('\n'):
    if Ei.strip():
        Ei = float(Ei.strip())
        ens_ryd.append(Ei)
ens_ryd = np.array(ens_ryd)
ens_eV = np.multiply(ens_ryd, 13.605698)
esp_avg = np.average(ens_ryd)
esp_std = np.std(ens_ryd)
print(esp_avg)
print esp_std
outputing = E_total, np.array(E_total + ens_eV)
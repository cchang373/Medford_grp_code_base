# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 14:30:48 2017

@author: benjamin
"""
import os
from pickle import dump

#import matplotlib
#import matplotlib.pyplot as plt
import numpy as np

#matplotlib.use('Agg')


def get_BEEF_ensemble(log):
#    cur_dir = os.getcwd()
#    print cur_dir
#    os.chdir(os.path.join(*log.split(os.sep)[0:-2]))
#    if os.path.isfile('../converged.log'):
        f = open(log)
        txt = f.read()
        f.close()
        _,E_total = txt.rsplit('total energy              =',1)
        E_total,_ = E_total.split('Ry',1)
        E_total = float(E_total.strip())
        E_total *= 13.605698
        try:
            _, ens = txt.rsplit('BEEFens 2000 ensemble energies',1)
            ens,_ = ens.split('BEEF-vdW xc energy contributions',1)
            ens.strip()
            ens_ryd = []
            for Ei in ens.split('\n'):
                if Ei.strip():
                    Ei = float(Ei.strip())
                    ens_ryd.append(Ei)
            ens_ryd = np.array(ens_ryd)
            ens_eV = np.multiply(ens_ryd, 13.605698)
            print(log)
#	    os.chkdir(cur_dir)
            return [E_total, np.array(E_total + ens_eV)]
        except:
            print(log)
#	    os.chdir(cur_dir)
            return E_total
    #        return E_total, np.array(np.add(E_total,ens_eV))


#get_BEEF_ensemble('/nv/hp13/bcomer3/shared/espresso_rutile/espresso_Rutile/molecules/CO/esp.log/log')

ensembles = {}
for root, dirs, files in os.walk(".", topdown=False):
    for name in files:
	esp_pres = [x for x in root.split(os.sep) if x=="esp.log"]
	non_vib = esp_pres==["esp.log"]
        if name == "log" and non_vib==True:
#	    print root.split(os.sep)
            ens_name = root.split(os.sep)[-2]+'_'+root.split(os.sep)[-3]#+'_'+root.split(os.sep)[-4]
            ensembles[ens_name] = get_BEEF_ensemble(os.path.join(root, name))
os.chdir('/nv/hp13/bcomer3/shared/espresso_rutile/espresso_Rutile/pickle_barrel')
#print(ensembles['H2_molecules'])
for key in ensembles:
    f = open(key,'w+')
#    print key
#    print ensembles[key]
    dump(ensembles[key],f)
    f.close()
#os.chdir(cur_dir)

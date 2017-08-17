# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 14:30:48 2017

@author: benjamin
"""
import os
from pickle import dump
import sys
sys.path.insert(0,"/home/benjamin/catmap")
sys.path.insert(0,"/home/benjamin/espresso/espresso_rutile/tools")
sys.path.insert(0,"/gpfs/pace1/project/chbe-medford/medford-share/users/bcomer3/espresso_rutile/tools")
from ThemoChem import surface_vibrations, gas_vibrations
import numpy as np



def get_BEEF_ensemble(log):
     if '2x3' in log:
         return 0
     f = open(log)
     txt = f.read()
     f.close()
     print(log)
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
#	    os.chkdir(cur_dir)
            return [E_total, np.array(E_total + ens_eV)]
     except:
            print(log)
            return E_total



ensembles = {}
cur_dir = os.path.abspath('.')
for root, dirs, files in os.walk(".", topdown=False):
    for name in files:
        esp_pres = [x for x in root.split(os.sep) if x=="esp.log"]
        non_vib = esp_pres==["esp.log"]
        slab_conv_pres = [x for x in root.split(os.sep) if x=="slab_convergence"]
        non_save = [x for x in root.split(os.sep) if x=="calc.save"] == []
        non_slab_conv = slab_conv_pres ==[]
        if name == "log" and non_vib==True and non_slab_conv == True:
#	    print root.split(os.sep)
            DFT_references = [x for x in root.split(os.sep) if x=='DFT_references']
            B3LYP_garb = [x for x in root.split(os.sep) if x=='B3LYP_molecules']
            if not DFT_references:
                ens_name = root.split(os.sep)[-3]+'-'+root.split(os.sep)[-2]+'-'+root.split(os.sep)[-4]
                ensembles[ens_name] = get_BEEF_ensemble(os.path.join(root, name))
            elif not B3LYP_garb:
                ens_name = root.split(os.sep)[-2]+'-'+root.split(os.sep)[-3]+'-'+root.split(os.sep)[-4]
                ensembles[ens_name] = get_BEEF_ensemble(os.path.join(root, name))
try:
    os.chdir('/nv/hp13/bcomer3/shared/espresso_rutile/espresso_Rutile/pickle_barrel')
except:
    try:
            os.chdir('/home/benjamin/espresso/espresso_rutile/espresso_Rutile/pickle_barrel')
    except:
        print('cannot find the pickles')

for key in ensembles:
    f = open(key,'w+')
    dump(ensembles[key],f)
    f.close()

os.chdir(cur_dir)

vibrations = {}
cur_dir = os.path.abspath('.')
for root, dirs, files in os.walk(".", topdown=False):
    for name in files:
        esp_pres = [x for x in root.split(os.sep) if x=="esp.log"]
        non_vib = esp_pres==["esp.log"]
        slab_conv_pres = [x for x in root.split(os.sep) if x=="slab_convergence"]
        non_save = [x for x in root.split(os.sep) if x=="calc.save"] == []
        non_slab_conv = slab_conv_pres ==[]
        if name == "vib.txt" and non_vib==False and non_slab_conv == True:
            print os.path.join(root,name)
            DFT_references = [x for x in root.split(os.sep) if x=='DFT_references']
            if not DFT_references:
                ens_name = root.split(os.sep)[-3]+'-'+root.split(os.sep)[-2]+'-'+root.split(os.sep)[-4]
                vibrations[ens_name] = surface_vibrations(os.path.join(root, name),300)
            else:
                ens_name = root.split(os.sep)[-2]+'-'+root.split(os.sep)[-3]+'-'+root.split(os.sep)[-4]
                vibrations[ens_name] = gas_vibrations(os.path.join(root, name),300,101325)



try:
    os.chdir('/nv/hp13/bcomer3/shared/espresso_rutile/espresso_Rutile/vib_pickle_barrel')
except:
    try:
            os.chdir('/home/benjamin/espresso/espresso_rutile/espresso_Rutile/vib_pickle_barrel')
    except:
        print('cannot find the vibrating pickles')
for key in vibrations:
    f = open('vib-'+key,'w+')
    dump(vibrations[key],f)
    f.close()

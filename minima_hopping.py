# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 12:20:20 2017

@author: benjamin
"""

from ase import io
import numpy as np
from ase.constraints import Hookean
from espresso import espresso

from ase.optimize.minimahopping import MinimaHopping
bond_dict = {
             ('N','O'): 1.44,
             ('H','O'): 0.96,
             ('H','N'): 1.01,
             ('N','N'): 1.46,
             ('O','O'): 1.48,
             ('H','H'): 0.74
            


}

def is_nonmetal(atom):
    if atom.symbol in ['N','O','H']:
        return True
    else:
        return False
def is_bond(atom1,atom2,input_return=False):
    d = atom1.position-atom2.position
    d = np.linalg.norm(d)
    sym1 = atom1.symbol
    sym2 = atom2.symbol
    dict_inp = tuple(sorted([sym1,sym2]))
    ro = bond_dict[dict_inp]
    if (ro-0.5)<d<(ro+0.2):
        if input_return==True:
            return ro
        else:
            return True
    else:
        return False
    

atoms = io.read('converged_slab.traj')
cons = atoms.constraints
for i in atoms:
    for j in atoms:
        if is_nonmetal(i) and is_nonmetal(j):
            if i.symbol==j.symbol!='O'or i.symbol!=j.symbol:
                if j.index>i.index:
                    if is_bond(i,j):
                        print(i.symbol,j.symbol)
                        ro = is_bond(i,j,input_return=True)
                        cons.append(Hookean(a1=i.index,a2=j.index,rt=ro,k=5))
calcargs = dict(xc='BEEF-vdW',
        kpts=(4, 4, 1), #only need 1 kpt in z-direction
        pw=400.,
        dw=4000.,
        spinpol=True,
        beefensemble=True,
        printensemble=True,
        convergence={'energy':1e-6,
                    'mixing':0.05,
                    'maxsteps':1000,
                    'diag':'david'},
        startingwfc='atomic',
        smearing='fd', #fermi-dirac electron smearing
        sigma=0.1, #smearing width
        dipole={'status':True}, #dipole corrections True turns them on
        #parflags='-nk 2',
        outdir ='esp.log')

calc= espresso(**calcargs)
atoms.set_calculator(calc)
init_mag = np.zeros(len(atoms))
init_mag[-1] = 0.1
atoms.set_initial_magnetic_moments(magmoms=init_mag)
atoms.set_constraints(cons)
opt = MinimaHopping(atoms=atoms,trajectory='hop.traj',restart='hop.pckl')
opt()

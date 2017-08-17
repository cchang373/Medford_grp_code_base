# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 17:01:10 2017

@author: benjamin
"""

from ase.optimize.minimahopping import MinimaHopping
from ase.constraints import FixAtoms, Hookean
from ase import io
from espresso import espresso

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

calc = espresso(**calcargs)

atoms = io.read('converged_slab.traj')
cons = atoms.constraints
atoms.set_calculator(calc)
init_mag = np.zeros(len(atoms))
init_mag[-1] = 0.1
atoms.set_initial_magnetic_moments(magmoms=init_mag)
hop = MinimaHopping(atoms,
                    Ediff0=2.5,
                    T0=4000.)
hop(totalsteps=10)
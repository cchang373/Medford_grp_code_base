# -*- coding: utf-8 -*-
"""
Created on Thu May 11 16:29:45 2017

@author: benjamin
"""
from ase.phasediagram import Pourbaix, solvated
#from ase.phasediagram import Pourbaix
import numpy as np


refs = solvated('Ti')
energy_list =[
                ('Ag',0/4),
                ('AgO',2.15/4),
                ('AgO',4.46/6),
                ('AgO',5.39/4),
                ('AgOH',3.28/6),
                ('AgOH',0.93/6),
                ('AgOH',2.03/6),
]
refs+=energy_list
#refs += [('Zn', 0.0), ('ZnO', -3.323), ('ZnO2(aq)', -2.921)]
print refs
pb = Pourbaix(refs, Ag=1, O=2)
U = np.linspace(-2, 2, 200)
pH = np.linspace(-2, 16, 300)
d, names, text = pb.diagram(U, pH, plot=True)

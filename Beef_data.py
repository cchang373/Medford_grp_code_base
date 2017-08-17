# -*- coding: utf-8 -*-
"""
Written by Ben Comer 8/9/2017
ben.comer@gatech.edu/benmcomer@gmail.com
"""

import numpy as np
from copy import deepcopy

class beef:
    """
    this is an energy data type, you must feed it an iterable variable of 
    with the first object being a floating point number and the second being a
    numpy array. It has all the normal operations such as adding and 
    subtracting, with functionality to return averages and standard devations
    quickly.
    +=, -=, *=, and /= are implemented
    
    adding and subtracting a constant adds and subtracts both the energy and 
    beef ensembles. adding or subtracting a beef ensemble adds/subtracts from
    the individual components the way you'd like it to.
    
    kwargs:
    
    energy_ens_list: a beef energy and beef ensemble in the form of a two item
    list, the list should have the beef energy in the first position and
    beef ensemble in the second e.g. [-48.39,np.array(-48.29,-47.89,...)]
    
    vib: can store vibrational frequencies, empty by default, must be passed
    in as a list
    
    path: allows you to attach a filepath to 
    
    
    
    
    
    """
    def __init__(self,energy_ens_list):
        eng_is_float = type(energy_ens_list[0])==float or type(energy_ens_list[0])==int or type(energy_ens_list[0]) ==np.float64
        ens_is_np_array = type(energy_ens_list[1])==np.ndarray
        if eng_is_float and ens_is_np_array:
            self.energy = energy_ens_list[0]
            self.ensemble = energy_ens_list[1]
        else:
            raise TypeError('make sure the you\'re feeding beef a list containing a float/int and a numpy array')
    #functions to pull stuff out
    def eng(self):
        return self.energy
    def ens(self):
        return self.ensemble
    def ensemble_avg(self):
        avg = np.mean(self.ensemble)
        return avg
    def ensemble_std(self):
        return np.std(self.ensemble)
    def error_bar(self):
        return self.ensemble_std()
    def list_form(self):
        return [self.energy,self.ensemble]
    def array_form(self):
        return [self.energy,self.ensemble]        
    def pull_vibs(self):
        return self.vib
    #This part sets up all the basic operations to make the beef object 
    #act like a float
    def __radd__(self,other):
        if isinstance(other,beef):
            return beef([self.energy+other.energy,self.ensemble-other.ensemble])
        else:
            return beef([self.energy+other,self.ensemble+other])
    def __add__(self, other):
        if isinstance(other,beef):
            return beef([self.energy+other.energy,self.ensemble+other.ensemble])
        else:
            return beef([self.energy+other,self.ensemble+other])
    def __sub__(self,other):
        if isinstance(other,beef):
            return beef([self.energy-other.energy,self.ensemble-other.ensemble])
        else:
            return beef([self.energy-other,self.ensemble-other])
    def __rsub__(self,other):
        if isinstance(other,beef):
            return beef([self.energy-other.energy,self.ensemble-other.ensemble])
        else:
            return beef([self.energy-other,self.ensemble-other])
    def __iadd__(self,other):
        if isinstance(other,beef):
            self.energy += other.energy
            self.ensemble +=other.ensemble
        else:
            self.energy += other
            self.ensemble += other
        return self
    def __isub__(self,other):
        if isinstance(other,beef):
            self.energy -= other.energy
            self.ensemble -=other.ensemble
        else:
            self.energy -= other
            self.ensemble -= other
        return self
    def __mul__(self,other):
        if isinstance(other,beef):
            return beef([self.energy*other.energy,self.ensemble*other.ensemble])
        else:
            return beef([self.energy*other,self.ensemble*other])
    def __truediv__(self,other):
        return beef([self.energy/other,self.ensemble/other])
    def __imul__(self,other):
        if isinstance(other,beef):
            self.energy *= other.energy
            self.ensemble *=other.ensemble
        else:
            self.energy *= other
            self.ensemble *= other
    def __itruediv__(self,other):
        if isinstance(other,beef):
            self.energy /= other.energy
            self.ensemble /=other.ensemble
        else:
            self.energy /= other
            self.ensemble /= other
    def __pow__(self,other):
        if isinstance(other,beef):
            print('this is not implemented because there is no reason to do it')
        else:
            return beef([self.energy**other,self.ensemble**other])
    #these make the object pretend to be a list of length 2
    def __getitem__(self, indices):
        """
        This makes beef objects indexable as if they were [float,np.array()]
        
        so if b is a beef object, this makes it so that b[0] returns a float 
        and b[1] returns the ensemble as a np.array()
        """
        return self.list_form()[indices]
    def __setitem__(self,key,value):
        if key == 0:
            self.energy = float(value)
        elif key == 1:
            if type(value)==np.ndarray:
                self.ensemble = value
            elif type(value)==list:
                self.ensemble = np.array(value)
            else:
                raise TypeError('ensemble must be list or numpy array')
    def __len__(self):
        """
        the length is always 2 since it's pretending to be [float,np.array()]
        """
        return 2
    def s(self):
        """
        prints out a summary of the object, the Beef energy, the ensemble 
        average, and the ensemble standard deviation in tuple form
        """
        print((self.energy,self.ensemble_avg(),self.ensemble_std()))
    def c(self):
        """
        returns a deep copy of the beef object
        """
        return beef(deepcopy(self.list_form()))
    def __copy__(self):
        """
        returns a deepcopy of the beef object
        """
        return beef(deepcopy(self.list_form()))
    #Functions to store infromation about the calculation
    def add_vibrations(self,vibs,overwrite=False):
        
        if type(vibs) ==list:
            if self.vib==[] and not overwrite:
                self.vib = vibs
            elif overwrite:
                self.vib = vibs
            else:
                print('object already had vibrational frequencies and overwrite us False, using previous frequencies')
        else:
            raise TypeError('vibrational frequencies must be in a list object')
    def add_path(self,log):
        self.path = log
    def clean_vibs(self,Cutoff=50): #cutoff just arbitrarily set to 50 cm^-1
        clean_vibs = []        
        for item in self.vib:
            try:
                l = float(item)
                clean_vibs.append(l)
            except:
                if type(item)==complex:
                    clean_vibs.append(np.real(item))
                if type(item) == str:
                    clean_vibs.append(50)#if theres some string just make it 50
        self.vib=clean_vibs
    
            
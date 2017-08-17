# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 08:50:27 2017

@author: benjamin
"""
import os
from tempfile import mkstemp
from shutil import move
from os import remove, close

#from subprocess import call
cur_dir=os.path.relpath(".","..")
import os
textToSearch = 'spinpol=True,\n'

textToReplace = 'dw=4000.,\n        spinpol=True,\n        beefensemble=True,\n        printensemble=True,\n'

def replace(file_path, pattern, subst):
    #Create temp file
    fh, abs_path = mkstemp()
    with open(abs_path,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))
    close(fh)
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)

for root, dirs, files in os.walk(".", topdown=False):
    for name in files:
        if name == "qn_opt.py":
            f = open(os.path.join(root, name))
            txt = f.read()
            f.close()
            ens_in = True
            try:
                _,__ = txt.rsplit('beefensemble=True',1)
            except:
                ens_in = False

            if ens_in == False:
#               tempFile = open( os.path.join(root, name), 'r+' )
                print(os.path.join(root, name))
                replace(os.path.join(root, name),textToSearch, textToReplace)
#               tempFile.close()
        else:
            pass
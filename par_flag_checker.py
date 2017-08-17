# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 09:23:57 2017

@author: benjamin
"""
import os
from tempfile import mkstemp
from shutil import move
from os import remove, close

cur_dir=os.path.realpath(".")

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
        if name == "qn_opt.py" or name == "":
            os.chdir(os.path.join(root))
            if os.path.isfile('run.sh') and os.path.isfile('qn_opt.py'):
#		print(os.path.join(root))
                f = open('run.sh')
                txt1 = f.read()
                f.close()
                g = open('qn_opt.py')
                txt2 = g.read()
                g.close()
                _,run_cores = txt1.rsplit('#PBS -l nodes=',1)
                run_cores,_ = run_cores.split(':ppn=',1)
                run_cores = run_cores.strip()
#                print(run_cores)
                _,opt_cores = txt2.rsplit('        parflags=\'-nk ',1)
                opt_cores,_ = opt_cores.split('\',\n        outdir =\'esp.log\')',1)
                opt_cores = opt_cores.strip()
#                print(opt_cores)
                if run_cores == opt_cores:
                    pass
                else:
                    print(os.path.join(root, name))
                    replace(os.path.join(name),'\'-nk '+opt_cores+'\'','\'-nk '+run_cores+'\'')
            os.chdir(cur_dir)

for root, dirs, files in os.walk(".", topdown=False):
     for name in files:
        if name == "vib_calc.py":
#	    print os.path.join(root)
            os.chdir(os.path.join(root))
            if os.path.isfile('run_vib.sh') and os.path.isfile('vib_calc.py'):
#                print(os.path.join(root))
                f = open('run_vib.sh')
                txt1 = f.read()
                f.close()
                g = open('vib_calc.py')
                txt2 = g.read()
                g.close()
                _,run_cores = txt1.rsplit('#PBS -l nodes=',1)
                run_cores,_ = run_cores.split(':ppn=',1)
                run_cores = run_cores.strip()
#                print(run_cores)
                _,opt_cores = txt2.rsplit('        parflags=\'-nk ',1)
                opt_cores,_ = opt_cores.split('\',\n        outdir =\'esp.log\')',1)
                opt_cores = opt_cores.strip()
#                print(opt_cores)
                if run_cores == opt_cores:
                    pass
                else:
                    print(os.path.join(root, name))
                    replace(os.path.join(name),'\'-nk '+opt_cores+'\'','\'-nk '+run_cores+'\'')
            os.chdir(cur_dir)


# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 09:42:27 2017

@author: benjamin
"""
import os
import sys
sys.path.insert(0,"/nv/hp13/bcomer3/shared/espresso_rutile/tools")
from Change_Cores import custom_replace

def If_queue_allows(dir_to_check):
    save = [x for x in dir_to_check.split(os.sep) if x == "calc.save"]
    if not save ==[]:
        return False
    os.system("qstat -u bcomer3 >> qstat.txt")
    queue_pull = open('qstat.txt','r')
    queue_pull = queue_pull.read()
    #print(queue_pull)
    os.remove('qstat.txt')
    num_jobs = queue_pull.count('bcomer3')
    #print(num_jobs)
    os.system("qstat -u bcomer3 -f >> qstat.txt")
    queue_pull = open('qstat.txt','r')
    queue_pull = queue_pull.read()
    core_pull = queue_pull
    #print(queue_pull)
    os.remove('qstat.txt')
    queue_pull = queue_pull.split("PBS_O_WORKDIR=")
    core_pull = core_pull.split("Resource_List.procs = ")
    #print(len(queue_pull))
    work_dirs = []
    tot_cores = 0
    for i in range(1,num_jobs+1):
	#print core_pull[i]
        work_dir,_ = queue_pull[i].rsplit("PBS_O_HOST=")
        work_dir = work_dir.replace('\n','')
        work_dir = work_dir.replace('\t','')
        work_dir = work_dir.replace(',','')
        #print(work_dir)
        work_dirs.append(work_dir)
        core_raw,_ = core_pull[i].rsplit("Resource_List.walltime")
        cores = core_raw.split("procs =")
	#print cores
        tot_cores = tot_cores + float(cores[0])
    if dir_to_check in work_dirs:
        return False
    if tot_cores >= 200:
        try:
            custom_replace(dir_to_check+'/run.sh','walltime=48:00:00','walltime=24:00:00')
        except:
            try:
                custom_replace(dir_to_check+'/run_vib.sh','walltime=48:00:00','walltime=24:00:00')
            except:
                pass
        return True
    if tot_cores < 200:
        try:
            custom_replace(dir_to_check+'/run.sh','walltime=24:00:00','walltime=48:00:00')
        except:
            try:
                custom_replace(dir_to_check+'/run_vib.sh','walltime=24:00:00','walltime=48:00:00')
            except:
                pass
        return True

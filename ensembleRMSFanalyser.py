#Mahakaran Sandhu, RSC, 2019#

import matplotlib.pyplot as plt
import os
import numpy as np
import subprocess
import random 
import multiprocessing





############################################################################################REPLICATE RUNNING MODULES##############################################################################################

def replicates_setup(parentdir, p_opt, w_opt, t_opt, modelpath, datapath, cifpath, ligandNames, nreps=10):
    """Ligand names is a list of names of ligands to restrain."""

    replicate_dirs = list(range(1,nreps+1))

    restr_sel = ['resname '+i +' or' if ligandNames.index(i)!=ligandNames.index(ligandNames[-1]) else 'resname '+i for i in ligandNames]
    restr_str = ' '.join(restr_sel)
    
    



    for i in replicate_dirs:
        random_seed = random.randint(1000000, 9000000)
        os.chdir(parentdir)
        dirname = 'replicate_' + str(i)
        outfilename = 'refine_' + str(i) +  '.sh'
        os.mkdir(dirname)
        filepath = os.path.join(dirname,outfilename)
        outfile = open(filepath, 'w+')
        outfile.write('#!/bin/bash \n \n')
        outfile.write('p=' + str(p_opt) + '\n')
        outfile.write('w=' + str(w_opt) + '\n')
        outfile.write('t=' + str(t_opt) + '\n \n')
        outfile.write('model='+modelpath+' \n')
        outfile.write('data='+datapath+' \n')
        outfile.write('cif='+cifpath+' \n \n')
        outfile.write("phenix.ensemble_refinement $model $data $cif harmonic_restraints.selections='"+str(restr_str)+"' ptls=$p wxray_coupled_tbath_offset=$w tx=$t  nproc=1 seed=" + str(random_seed))
        outfile.close()



 
def runEnsembleRefinement(dirname):
    print ('Now running:'+str(dirname))
    os.system('cd '+ str(dirname) + '&& sh *sh')
    print ('Finished running:'+str(dirname))

def chunkIt(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out

def run_replicates(parentdir, cores):
    """For a list of directories, use cores number of cores to run the ER protocol"""
    replicate_dirs = [i for i in os.listdir(parentdir) if os.path.isdir(i) and 'replicate' in i]
    completedER = []
    
    
    x = len(replicate_dirs)/cores
    
    
    
    processedDirList = chunkIt(replicate_dirs, x)
    print (processedDirList)

    for i in processedDirList:
        procs = []
        for j in i:
            print(j)
            proc = multiprocessing.Process(target=runEnsembleRefinement, args=(j,))
            procs.append(proc)
            proc.start()
        for k in procs:
            k.join()









#################################################################################################################ANALYSIS MODULES##############################################################################################################


def split_chains(list_of_tuples):
    resids = [int(list_of_tuples[i][0]) for i in range(len(list_of_tuples))]
    breaks = []
    for i in range(len(resids)):
        if (len(str(resids[i-1]))!=len(str(resids[i]))):
            if (resids[i]!=100) and i!=0 :
                breaks.append(i)
    breaks.append(-1)
    
    chains = []
    for i in range(len(breaks)):
        if i==0:
            chain = list_of_tuples[:breaks[i]]

        else:
            chain = list_of_tuples[breaks[i-1]:breaks[i]]
        chains.append(chain)
    
    return chains
        



def processrmsf_1(xvg_file):
    """Takes an input xvg file (string) with rmsfs"""

    data = open(xvg_file)
    data = data.readlines()
    data = [i.strip() for i in data]
    data = [i for i in data if ('@' not in i) and ('#' not in i)]
    proc_data = []

    for i in data: 
        r = i.split()
        proc_data.append(r)
    
    chain_data = split_chains(proc_data)
 
    return chain_data




def processrmsf_2(inlist):

    data_series=[]
    for i in inlist:
        t_array = np.array(i)
        x_data=(list(t_array[:,0]))
        y_data=(list(t_array[:,1]))
        x_data=[int(i) for i in x_data]
        y_data=[float(i) for i in y_data]

        series = [x_data,y_data]

        data_series.append(series)
        
    return data_series


def processrmsf_3(inlist):
    chains_resids = []
    for i in inlist:
        chains_resids.append(i[0])

    average_resids = list(set.intersection(*map(set,chains_resids)))

   
    to_average_data = []
    for resid in average_resids:
        res_data_to_average = [resid]
        for chain in inlist:
            data_index = chain[0].index(resid)
            rmsf_dat = chain[1][data_index]
            res_data_to_average.append(rmsf_dat)
        to_average_data.append(res_data_to_average)

    return to_average_data            



def processrmsf_4(inlist):
    averaged_data = []
    for i in inlist:
        res_id = i[0]
        rmsf_data = i[1:]
        av_rmsf = np.mean(rmsf_data)
        av_list = [res_id,av_rmsf]
        #print (av_list)
        averaged_data.append(av_list)
    return averaged_data


def plot_xvg(xvg_file):
    """Takes an input xvg file (string), label for the x-axis (string) and for the y-axis (string);
    plots a plot"""

    data = open(xvg_file)
    data = data.readlines()
    data = [i.strip() for i in data]
    data = [i for i in data if ('@' not in i) and ('#' not in i)]
    x_axis=[]
    y_axis = []
    for i in data: 
        r = i.split()
        x_axis.append(float(r[0]))
        y_axis.append(float(r[1]))
        
    return (x_axis, y_axis)
    




def plot_replicate_RMSFs(parentdir):

    ER_directory = parentdir
    subplt_dims = (5,2)
    figure,axs = plt.subplots(subplt_dims[0], subplt_dims[1], figsize=(15,15))

    rows = [i for i in range(subplt_dims[0])]
    cols = [i for i in range(subplt_dims[1])]

    sbplt_size = len(rows)*len(cols)

    j=0
    k=0
    for i in range(1,sbplt_size+1):
        try: 
            row_num = rows[j]
            col_num = cols[k]
        except:
            j=0
            k+=1
            row_num = rows[j]
            col_num = cols[k]
        j+=1    

        try:
            rmsf_path = ER_directory+'/replicate_'+str(i)+'/rmsf.xvg'
            rmsf_data = plot_xvg(rmsf_path)
            axs[row_num,col_num].plot(rmsf_data[0], rmsf_data[1])
            axs[row_num,col_num].set_title('replicate_'+str(i))
        except:
            pass


        



def ensemble_RMSF(replicate_dir):
    return np.array(processrmsf_4(processrmsf_3(processrmsf_2(processrmsf_1(replicate_dir+'rmsf.xvg')))))






def find_min_Rfree(parentdir):
    import glob
    os.chdir(parentdir)
    data_list = []
    rep_dirs = [i for i in os.listdir() if os.path.isdir(i) and ('replicate' in i) ]
    for i in rep_dirs:
        os.chdir(str(i))
        logfile = glob.glob("*.log")
        logfile = open(str(logfile[0]), 'r')
        lines = logfile.readlines()
        for j in lines: 
            if ('FINAL' in j) and ('Rfree' in j):
                r_free = float(j[28:36])
        rep_name = i
        out_tuple = (i,r_free)
        data_list.append(out_tuple)
        os.chdir(parentdir)
        r_frees = [data_list[i][1] for i in range(len(data_list))]
        repl_name = [data_list[i][0] for i in range(len(data_list))]
    print (repl_name)
    print (r_frees)
    return (data_list[r_frees.index(min(r_frees))][0], data_list[r_frees.index(min(r_frees))][1])




def get_rmsfs(parent_dir): 

    rep_dirs = [i for i in os.listdir() if os.path.isdir(i) and ('replicate' in i) ]
    
    for i in rep_dirs:
        os.chdir(i)
        print(i)
        subprocess.call('gunzip *.gz', shell=True)
        subprocess.call('grep -vwE "(HOH|MPD|CAC|ZN)" *ensemble.pdb > clean_ensemble.pdb', shell=True)
        subprocess.call("echo '3'|gmx rmsf -f clean_ensemble.pdb -s clean_ensemble.pdb -o rmsf.xvg -res", shell=True)
        os.chdir(parent_dir)

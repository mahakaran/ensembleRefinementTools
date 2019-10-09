import matplotlib.pyplot as plt
import os
import numpy as np



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
            print(data_index)
            rmsf_dat = chain[1][data_index]
            print(rmsf_dat)
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



def ensemble_RMSF(xvgfile):
    return processrmsf_4(processrmsf_3(processrmsf_2(processrmsf_1(xvgfile))))





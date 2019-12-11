#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 11:01:48 2017

@author: spiros
"""

import pickle
import os
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from place_cell_metrics import sparsity_index2
from place_cell_metrics import selectivity_index
from place_cell_metrics import peak_frequency
from place_cell_metrics import field_size
from tqdm import tqdm

fnames = 'Simulation_Results/'

Npyramidals = 130
Nperms      = 500 
Nbins       = 100

npath_x,npath_y = 200, 1
  
nTrials = 10


spec='final_results'

os.system('mkdir -p '+spec+'/figures_final/')

path_figs = spec+'/figures_final/'
path_data = spec+'/data_final/'

file_load = spec+'/metrics2/'
file_load_perms = spec+'/metrics_permutations/'
trials = [str(i) for i in range(1,nTrials+1)]
maindir=os.getcwd()

my_list = ['Control']
for i in range(5,36,5):
    my_list+= ['Desynch'+str(i), 'ALL0.50_'+str(i)]
for i in ['0.75','0.50','0.25']:
    my_list+= ['SOMred'+i, 'PVred'+i]
my_list+=['SOMdel', 'PVdel']

everything = {}

rateMaps = {}
infors   = {}
stabs1   = {}
stabs2   = {}
stabsAll = {}

prog_bar1 = tqdm(my_list)
prog_bar2 = tqdm(trials)

for case in prog_bar1:
    
    prog_bar1.set_description('Processing Case: %s' % case)
    
    infor      = []
    stab1      = []
    stab2      = []
    staball    = []
    for ntrial in prog_bar2:
        prog_bar2.set_description('Processing Trial %s' % ntrial)
        with open(file_load+'/pickled_sn_'+case+'_'+ntrial+'.pkl', 'rb') as f:
            loaded_data=pickle.load(f)    
                
        if ntrial=='1':
            rateMaps[case] = loaded_data['maps']
        else:
            rateMaps[case] = np.concatenate((rateMaps[case], loaded_data['maps']), axis=0)

        for npyr in xrange(Npyramidals):
            
            with open(file_load_perms+'perms_pickled_info_'+case+'_Npyr_'+str(npyr)+'_Mouse_'+ntrial+'.pkl', 'rb') as f:
                infor.append(pickle.load(f))
            with open(file_load_perms+'perms_pickled_stab_even_odd_'+case+'_Npyr_'+str(npyr)+'_Mouse_'+ntrial+'.pkl', 'rb') as f:
                stab1.append(pickle.load(f))
            with open(file_load_perms+'perms_pickled_stab_first_second_'+case+'_Npyr_'+str(npyr)+'_Mouse_'+ntrial+'.pkl', 'rb') as f:
                stab2.append(pickle.load(f))  
            with open(file_load_perms+'perms_pickled_stab_all_'+case+'_Npyr_'+str(npyr)+'_Mouse_'+ntrial+'.pkl', 'rb') as f:
                stab_all=pickle.load(f)  
            staball.append([np.nanmean(i) for i in stab_all])
        
    infors[case] = infor
    stabs1[case] = stab1
    stabs2[case] = stab2
    stabsAll[case] = staball
    
nCells = rateMaps[my_list[0]].shape[0]      

SpatialMapsALL =  {}
SpatialMapsPLC =  {} 
    
for case in prog_bar1:
    prog_bar1.set_description('Processing Case: %s' % case)

    inforALL1           = []
    stabALL1            = []
    inforALL2           = []
    stabALL2            = []
    sparsALL            = []
    selecALL            = []
    peaksALL            = []
    sizesALL            = []
    average1ALL         = []
    average2ALL         = []
    place_cells_idx     = []
    
    numbersALL          = []
    numbers_plc  = 0
    for npyr in xrange(nCells):
        
        pval_infor   = sum(infors[case][npyr][0]>=infors[case][npyr][1:])/float(Nperms)
        pval_stab1   = sum(stabs1[case][npyr][0]>=stabs1[case][npyr][1:])/float(Nperms)
        pval_stab2   = sum(stabs2[case][npyr][0]>=stabs2[case][npyr][1:])/float(Nperms)
        pval_staball = sum(stabsAll[case][npyr][0]>=stabsAll[case][npyr][1:])/float(Nperms)

        rate_map = rateMaps[case][npyr,:,:]
       
        maxpeak   = np.max(rate_map)
        meanrate  = np.mean(rate_map)
        fsize     = field_size(rate_map, relfreq=0.2*maxpeak, track_length=Nbins)[0]
        infield   = field_size(rate_map, relfreq=0.2*maxpeak, track_length=Nbins)[1]
        outfield  = field_size(rate_map, relfreq=0.2*maxpeak, track_length=Nbins)[2]
       
        clevel=0.99
        if (maxpeak>=1.0) and (pval_infor>=clevel) and (pval_staball>=clevel) and (5./(npath_x/Nbins)<=fsize<=25./(npath_x/Nbins)):# and (npyr in place_cells):#<=50./(npath_x/Nbins):        

            numbers_plc += 1
            place_cells_idx.append(npyr)

            inforALL2.append(infors[case][npyr][0])
            #stabALL2.append((stabs1[case][npyr][0]+stabs2[case][npyr][0])/2.0)
            stabALL2.append(stabsAll[case][npyr][0])
            sparsALL.append(sparsity_index2(rate_map))
            selecALL.append(selectivity_index(rate_map))
            peaksALL.append(peak_frequency(rate_map))
            sizesALL.append(fsize)
            average1ALL.append(infield)
            average2ALL.append(outfield)
        
        inforALL1.append(infors[case][npyr][0])
        #stabALL1.append((stabs1[case][npyr][0]+stabs2[case][npyr][0])/2.0)
        stabALL1.append(stabsAll[case][npyr][0])    
        
    numbersALL.append(numbers_plc/(len(trials)*130.0))


    rate_maps_all = rateMaps[case][:,:,:]
    Ncells = rate_maps_all.shape[0]
    
    idx        = np.argmax(rate_maps_all.squeeze(), axis=1)
    new_idx    = np.lexsort((range(Ncells), idx))
    rtMaps     = rate_maps_all[new_idx,:,:].squeeze()
    
    Max        = np.max(rtMaps, axis=1).reshape(-1,1)
    for i in xrange(Max.shape[0]):
        if Max[i,0]==0:
           Max[i,0]=1e-12 
              
    rtMaps  = rtMaps / np.repeat(Max, Nbins, axis=1) 
    SpatialMapsALL[case] = rtMaps
    
    fig, axes = plt.subplots(nrows=1, ncols=2)
    
    im0 = axes[0].imshow(rtMaps.squeeze(), cmap="jet", aspect='equal')
    divider = make_axes_locatable(axes[0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im0, cax=cax)                    
    ax = plt.gca()

    rate_maps_plc = rateMaps[case][place_cells_idx,:,:]
    Ncells = rate_maps_plc.shape[0]
    
    idx        = np.argmax(rate_maps_plc.squeeze(), axis=1)
    new_idx    = np.lexsort((range(Ncells), idx))
    rtMaps     = rate_maps_plc[new_idx,:,:].squeeze()
    
    Max        = np.max(rtMaps, axis=1).reshape(-1,1)
    for i in xrange(Max.shape[0]):
        if Max[i,0]==0:
           Max[i,0]=1e-12 
              
    rtMaps  = rtMaps / np.repeat(Max, Nbins, axis=1) 
    
    if rtMaps.shape[0] != 0:
        im1 = axes[1].imshow(rtMaps.squeeze(), cmap="jet", aspect='equal')
        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        divider = make_axes_locatable(axes[1])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im1, cax=cax)
    
    SpatialMapsPLC[case] = rtMaps
    
    plt.tight_layout()
    plt.savefig(path_figs+'/'+case+'_PlaceCells_heatmap.pdf',format='pdf',dpi=300)
    plt.savefig(path_figs+'/'+case+'_PlaceCells_heatmap.png',format='png',dpi=300)
    plt.cla()
    plt.clf()
    plt.close()

    mydict_all={}                    
    mydict_all['information_all']    = inforALL1
    mydict_all['stabilityIndex_all'] = stabALL1
    mydict_all['information_plc']    = inforALL2
    mydict_all['stabilityIndex_plc'] = stabALL2
    mydict_all['sparsity']           = sparsALL
    mydict_all['selectivity']        = selecALL
    mydict_all['fieldsize']          = sizesALL
    mydict_all['PeakRate']           = peaksALL
    mydict_all['numbers_plc']        = numbersALL
    mydict_all['MeanRateInfield']    = average1ALL
    mydict_all['MeanRateOutfield']   = average2ALL

    metrics = ['information_all', 'stabilityIndex_all', 'information_plc', 'stabilityIndex_plc',
               'sparsity', 'selectivity', 'fieldsize', 'PeakRate', 'MeanRateInfield', 'MeanRateOutfield']

    everything[case] = mydict_all

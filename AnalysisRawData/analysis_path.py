#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 21:28:15 2017

@author: spiros
"""

import numpy as np
import pickle, sys, os, time
import scipy.ndimage.filters as flt
from functions_analysis import spike_map, binning
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def analysis_path_cluster(ntrial, case):

    folder1 = 'final_results'
    os.system('mkdir -p '+folder1+'/metrics/')
    fdname1 = '/'+folder1+'/figures2/'
    fdname2 = '/'+folder1+'/metrics2/'
    
    print "Analyse ... " + case + " Trial " +ntrial
    
    os.system('mkdir -p '+folder1+'/figures2/')
    os.system('mkdir -p '+folder1+'/metrics2/')
    maindir=os.getcwd()

    
    # Give path dimensions
    npath_x, npath_y = 200, 1
    # Number of pyramidal
    Ncells  = 130
    Nbins   = 100
    runsAll = 10
    
    # Gaussian filter parameters
    sigma_c    = 5.0/(npath_x/Nbins)
    truncate_c = 30.0/(npath_x/Nbins)
    
    ### Define the map size!
    # 3-d matrix of all pyramidals
    rateMaps = np.zeros((Ncells,Nbins,npath_y))
    rateMaps_unsmoothed = np.zeros((Ncells,Nbins,npath_y))
    time_array_in_bin = np.zeros((Ncells,Nbins,npath_y))
    
    # File location - pathfile
    fileload  = 'peyman_results_new'+inum +'/metrics_permutations/'
    
    with open(fileload+'path_all_trial_'+str(ntrial)+'_'+case+'.pkl', 'rb') as f:
        path_all=pickle.load(f)    
    
    with open(fileload+'spiketimes_all_trial_'+str(ntrial)+'_'+case+'.pkl', 'rb') as f:
        spiketimes_all=pickle.load(f)
    
    # Loop for all pyramidals
    for npyr in range(Ncells):
        # A matrix for rate map
        Zall = np.zeros((Nbins,npath_y))
        time_array_all = np.zeros(Nbins*npath_y)

        for nrun in range(1,runsAll+1):
            
            # Load of path -- different for each run
            path = path_all[nrun-1]

            time_array = np.bincount(path[:,0])[1:]
            csum = np.cumsum(time_array)

            spiketimes = spiketimes_all['Run'+str(nrun)]['Pyramidal'+str(npyr)][case]

            Z       = spike_map(spiketimes,csum,npath_x,npath_y)
            Zbinned = binning(Z, Nbins, 'summing')
            time_binned = binning(time_array, Nbins, 'summing').squeeze()
            
            # Take the sum over all runs given by total
            Zall += Zbinned
            time_array_all += time_binned # time spent in each bin in ms

        # Calculate the time spent in each bin
        # convert to Hz, so divide with seconds,  time ms/1000 (ms/sec) --> seconds
        time_array_sec = (time_array_all/1000.0).reshape(Zall.shape)
        
        # Gaussian smoothing
        time_array_fil = flt.gaussian_filter1d(time_array_sec, axis=0, sigma=sigma_c, 
                                               mode='nearest',truncate = truncate_c)
        
        Zsmoothed=flt.gaussian_filter1d(Zall, axis=0,sigma=sigma_c, 
                                        mode='nearest',truncate = truncate_c)
        
        
        Zmean = np.divide(Zsmoothed, time_array_fil)
        
        # Spatial Coherence - "Spatial representations of place cells in darkness are supported by path integration and border information"

        rateMaps_unsmoothed[int(npyr),:,:]=Zall
        rateMaps[int(npyr),:,:]=Zmean
        time_array_in_bin[int(npyr),:,:]=time_array_fil
        

    print '\nDone with the rate maps'

    fig, axes = plt.subplots(nrows=13, ncols=10,figsize=(15, 15))
    nn=0
    for ax in axes.flat:
        Max = np.max(rateMaps[nn,:,:])
        im = ax.imshow(rateMaps[nn,:,:].T/Max, origin='lower',cmap="jet", aspect='10')
        ax.tick_params( axis='y',which='both',right='off',left='off',labelleft='off')
        ax.title.set_text('PC'+str(nn) + '  ' + str(np.round(Max,1))+ ' Hz')
        nn+=1

    fig.colorbar(im, ax=axes.ravel().tolist())
    
    if not os.path.exists(maindir+fdname1+'/'):
        os.makedirs(maindir+fdname1+'/')    
        
    plt.savefig(maindir+fdname1+'/'+case+'_'+ntrial+'_heatmap.pdf',format='pdf',dpi=300)
    
    idx        = np.argmax(rateMaps.squeeze(), axis=1)
    new_idx    = np.lexsort((range(Ncells), idx))
    rtMaps     = rateMaps[new_idx,:,:].squeeze()
    
    Max        = np.max(rtMaps, axis=1).reshape(-1,1)
    for i in xrange(Max.shape[0]):
        if Max[i,0]==0:
           Max[i,0]=1e-12 
                
    fig = plt.subplots(figsize=(15, 15))
    ax = plt.gca()
    im = ax.imshow(rtMaps, cmap="jet", aspect='equal')
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    ax.set_xlim((0, Nbins))
    ax.set_xticks(range(0,Nbins+1,Nbins/4))
    ax.set_xticklabels([str(x) for x in range(0,npath_x+1,npath_x/4)], fontsize = 13)
    ax.set_yticks(range(0,Ncells+1,20))
    ax.set_yticklabels([str(x) for x in range(0,Ncells+1,20)], fontsize = 13)    
    ax.set_title(case, fontsize=14)
    plt.savefig(maindir+fdname1+'/'+case+'_'+ntrial+'_heatmap_all_cells.pdf',format='pdf',dpi=300)
    #==============================================================================
    # ##################### RATE MAPS SAVING #################################
    #==============================================================================

    mydict= {}
    mydict['maps']=rateMaps
    mydict['maps_unsmoothed']=rateMaps_unsmoothed
    mydict['time_in_bin'] = time_array_in_bin
    
    filesave = maindir+fdname2
    if not os.path.exists(filesave):
        os.makedirs(filesave)
    
    with open(filesave+'/pickled_sn_'+case+'_'+ntrial+'.pkl', 'wb') as handle:
        pickle.dump(mydict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    
    print "\nDone with "+case+" analysis. Done with trial "+ntrial


tic      = time.time()
ntrial   = sys.argv[1]
case     = sys.argv[2]
results  = analysis_path_cluster(ntrial, case)
toc      = time.time()

print "\nTotal time: "+str(round(toc-tic,3))+" seconds"

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 14 10:36:49 2018

@author: spiros
"""

def permutations_shuffle(case,npyr,path_all, spiketimes_all, npath_x, npath_y,Nperms, Nbins, nruns):
    
    import numpy as np
    import scipy.ndimage.filters as flt
    from functions_analysis import spike_map, binning
    from place_cell_metrics import spatial_info, stability_index

    # Gaussian filter parameters
    sigma_c    = 5.0/(npath_x/Nbins)
    truncate_c = 30.0/(npath_x/Nbins)
    
    spinfo                 = []
    stability_even_odd     = []
    stability_first_second = []
    stability_all          = []
    
#    np.random.seed(1)
    shift = np.random.randint(1, npath_x+1, (Nperms, nruns))
    shift[0,:] = 0
    for i in range(Nperms):
        print "Permutation: " + str(i)
        
        # A matrix for rate map
        Zall    = np.zeros((Nbins,npath_y))
        ZOdd    = np.zeros((Nbins,npath_y))
        ZEven   = np.zeros((Nbins,npath_y))
        Zfirst  = np.zeros((Nbins,npath_y))
        Zsecond = np.zeros((Nbins,npath_y))      
        Zrun    = np.zeros((nruns,Nbins,npath_y))
        
        time_array_all    = np.zeros(Nbins*npath_y)
        time_array_odd    = np.zeros(Nbins*npath_y)
        time_array_even   = np.zeros(Nbins*npath_y)
        time_array_first  = np.zeros(Nbins*npath_y)
        time_array_second = np.zeros(Nbins*npath_y)
        
        for nrun in range(1, nruns+1):
            
            ### Load PATH and SPIKETIMES of Pyramidal npyr
            # Load of path -- different for each run
            position = path_all[nrun-1][:,0]
            # Load Spiketimes for all PCs
            spiketimes = spiketimes_all['Run'+str(nrun)]['Pyramidal'+str(npyr)][case]
            
            # take a circular shuffle by a random number    
            # make circular shuffling - Daniel Aharoni
            pos  = shift[i, nrun-1]
            position_shift = position + pos
            position_new = np.mod(position_shift, npath_x)+1
            
            time_array_shifted = np.bincount(position_new)[1:]
            csum = np.cumsum(time_array_shifted)
            
            time_array_bin = binning(time_array_shifted, Nbins, 'summing').squeeze()
            # Convert to ms to seconds
            time_array_sec = time_array_bin / 1000.
            # Gaussian filter of the time signal
            time_array_sm = flt.gaussian_filter1d(time_array_sec, axis=0, sigma=sigma_c,
                                                  mode='nearest', truncate = truncate_c)
            
            Z = spike_map(spiketimes,csum,npath_x,npath_y)
                        
            Zbin = binning(Z, Nbins, 'summing')

            Zsm = flt.gaussian_filter1d(Zbin, axis=0, sigma=sigma_c, mode='nearest', truncate = truncate_c)
            Zrun[int(nrun-1),:,:] = np.divide(Zsm, (time_array_sm).reshape(Nbins,npath_y))
            
            
            # Odd - Even trials - runs
            if nrun % 2 != 0:
                ZOdd  = np.add(ZOdd, Zbin)
                time_array_odd = np.add(time_array_odd, time_array_sec)
            else:
                ZEven  = np.add(ZEven, Zbin)
                time_array_even = np.add(time_array_even, time_array_sec)
            # First/Second half trials - runs
            if nrun <= (nruns/2):
                Zfirst  = np.add(Zfirst, Zbin)
                time_array_first = np.add(time_array_first, time_array_sec)
            else:
                Zsecond = np.add(Zsecond, Zbin)
                time_array_second = np.add(time_array_second, time_array_sec)         
            
            # Take the sum over all runs given by total
            Zall = np.add(Zall, Zbin)
            time_array_all = np.add(time_array_all, time_array_sec) # time spent in each bin in seconds
    
        time_array_all_sm    = flt.gaussian_filter1d(time_array_all,axis=0, 
                                                  sigma=sigma_c, mode='nearest', truncate=truncate_c)
        time_array_odd_sm    = flt.gaussian_filter1d(time_array_odd, axis=0,
                                                  sigma=sigma_c, mode='nearest', truncate=truncate_c)
        time_array_even_sm   = flt.gaussian_filter1d(time_array_even, axis=0,
                                                  sigma=sigma_c, mode='nearest', truncate=truncate_c)
        time_array_first_sm  = flt.gaussian_filter1d(time_array_first, axis=0,
                                                  sigma=sigma_c, mode='nearest', truncate=truncate_c)
        time_array_second_sm = flt.gaussian_filter1d(time_array_second, axis=0,
                                                  sigma=sigma_c, mode='nearest', truncate=truncate_c)
        
        # Gaussian smoothing
        Zsmoothed     = flt.gaussian_filter1d(Zall, axis=0,
                                              sigma=sigma_c, mode='nearest', truncate=truncate_c)
        ZsmoothedOdd  = flt.gaussian_filter1d(ZOdd, axis=0,
                                              sigma=sigma_c, mode='nearest', truncate=truncate_c)
        ZsmoothedEven = flt.gaussian_filter1d(ZEven, axis=0,
                                              sigma=sigma_c, mode='nearest', truncate=truncate_c)
        Zsmoothed1    = flt.gaussian_filter1d(Zfirst, axis=0,
                                              sigma=sigma_c, mode='nearest', truncate=truncate_c)
        Zsmoothed2    = flt.gaussian_filter1d(Zsecond, axis=0,
                                              sigma=sigma_c, mode='nearest', truncate=truncate_c)
        
        # convert to Hz, so divide with seconds
        Zmean = np.divide(Zsmoothed, (time_array_all_sm).reshape(Nbins,npath_y))
        
        ZmeanOdd  = np.divide(ZsmoothedOdd, (time_array_odd_sm).reshape(Nbins,npath_y))
        ZmeanEven = np.divide(ZsmoothedEven, (time_array_even_sm).reshape(Nbins,npath_y))

        Zmean1  = np.divide(Zsmoothed1, (time_array_first_sm).reshape(Nbins,npath_y))
        Zmean2  = np.divide(Zsmoothed2, (time_array_second_sm).reshape(Nbins,npath_y))
             
        
        # calculate spatial information        
        spinfo.append(spatial_info(Zmean, time_array_all.reshape(-1,1)))
        
        # calculate stability
        stability_even_odd.append(stability_index(x=ZmeanOdd,y=ZmeanEven))
        stability_first_second.append(stability_index(x=Zmean1,y=Zmean2))
        stability_all.append(stability_index(x=Zrun, y=None))
        
        if i==0:
            rateMaps_all = Zmean
        
        mydict={}
        mydict['information']     = spinfo
        mydict['stabEvenOdd']     = stability_even_odd
        mydict['stabFirstSecond'] = stability_first_second
        mydict['stabAll']         = stability_all
        mydict['rateMaps']        = rateMaps_all
        
    return mydict

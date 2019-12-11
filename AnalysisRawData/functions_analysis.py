def path_analysis(path,npath_x,npath_y):
    import numpy as np
    '''
    A function that transforms the MATLAB's path to time vector
    and gives a full path of the subject.
    
    path is the loaded path.txt--output from matlab code (grid cells)
    '''
    # Path length
    L =path.shape[0]
    # Time spent at each point
    time_array = np.zeros(npath_x*npath_y)
    # Path -- for checking
    path_full = np.zeros((npath_x*npath_y,2))
    
    counter=1
    j=0
    z=0
    for i in range(1,L):
        if path[i,0]==path[i-1,0] and path[i,1] == path[i-1,1]:
            # counter is 1ms as MATLAB simulates with this time step
            counter+=1
            time_array[j]=counter
            # Create the path - each row is one point!
            path_full[z]=int(path[i,0]),int(path[i,1])
        else:
            j+=1
            counter=1 # 1ms in each entry
            z+=1
    
    # Cumulative sum of time array and add 0 to the beginning
    csum = np.cumsum(time_array)
    csum=np.insert(csum,0,0)
    
    return time_array, path_full, csum


def spike_map(spiketimes,csum, npath_x, npath_y):
    '''
    Make the Spike Matrix from spiketimes, path and time_array
    '''
    
    import numpy as np
    
    Z = np.zeros((npath_x,npath_y))
    
    if len(spiketimes) != 0:
        for spike in spiketimes:
            # Find the last positive index --> this is the mapping from time to space 
            
            if spike > csum[-1]:
                continue
            else:
                idxs = np.argwhere((spike - csum) > 0)
                if idxs.shape[0]==0:
                    idx=0
                else:
                    idx = idxs.shape[0]
                
    
                Z[idx,:] += 1
    
    if Z.shape[0] != npath_x or Z.shape[1] != npath_y:
        print 'Error in Z dimensions'
    
    return Z


#def spike_map(spiketimes,time_array,csum,npath_x,npath_y,path_full):
#    '''
#    Make the Spike Matrix from spiketimes, path and time_array
#    '''
#    
#    import numpy as np
#    
#    path_spiketimes = np.zeros((npath_x*npath_y))
#    for j in range(1,len(csum)):
#        counter=1
#        for spike in spiketimes:
#            if csum[j-1]<=spike<csum[j]:
#                path_spiketimes[j-1]=counter
#                counter+=1
#    
#    
#    # Loop over path_full to have correct
#    # path coordinates
#    Z = np.zeros((npath_x,npath_y))
#    idx=0
#    for ic in path_full:
#        xc = int(ic[0])-1
#        yc = int(ic[1])-1
#        Z[xc,yc]=path_spiketimes[idx] # number of spikes in each bin
#        idx+=1
#    
#    if Z.shape[0] !=npath_x or Z.shape[1] !=npath_y:
#        print 'Error in Z dimensions'
#    
#    return Z


def binning(a, N, method):
    import numpy as np
    
    if (isinstance(N, int)==False):
        raise ValueError('Binning size must be integer.')
        
    a = a.squeeze()
    L = a.shape[0]
    rem = L % N
    
    if rem !=0:
        raise ValueError('Not a valid binning size.')
    
    # find the step
    step = L/N
    
    b = np.zeros(N)
    
    cnt = 0    
    for i in range(0, L, step):
        if method=='summing':
            b[cnt] = np.sum(a[i:i+step])
        elif method=='mean':
            b[cnt] = np.mean(a[i:i+step])
        cnt+=1
    
    
    return b.reshape(-1,1)
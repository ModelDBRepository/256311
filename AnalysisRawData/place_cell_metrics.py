import numpy as np
import sys


def overall_mean(big_matrix, N):
    overallMean = 0
    for npyr in range(N):
        overallMean += np.mean(big_matrix[npyr, :, :])
    return overallMean


def spatial_info(rate_matrix, time_bin):
    if rate_matrix.shape != time_bin.shape:

        sys.exit('Error. Dimension mismatch')

    # Time per bin
    tall = np.sum(time_bin)

    # Mean firing rate
    pall = time_bin/tall
    l_ = np.matmul(pall.T, rate_matrix).item()
    if l_ == 0:
        l_ = 1e-15

    ratio = rate_matrix/l_

    ratio[ratio == 0] = 1e-15

    spatial_info = np.sum(np.multiply(
        pall, np.multiply(ratio, np.log2(ratio))))

    return spatial_info


def spatial_info_attila(rate_matrix, time_bin):

    # Time per bin
    tall = np.sum(time_bin)

    # Occupancy
    pall = time_bin/tall
    # Mean firing rate
    l_ = np.matmul(pall.T, rate_matrix).item()
    if l_ == 0:
        l_ = 1e-15

    ratio = rate_matrix/l_
    ratio[ratio == 0] = 1e-15

    spatial_info = np.sum(np.multiply(
        pall, np.multiply(rate_matrix, np.log2(ratio))))

    return spatial_info


def selectivity_index(rate_matrix):
    A = np.mean(rate_matrix)
    B = np.max(rate_matrix)
    if A == 0:
        A = 1e-15

    selectivity = float(B)/A

    return selectivity


def sparsity_index(rate_matrix, time_bin, probability='uniform'):
    # Size of field
    n1 = rate_matrix.shape[0]
    n2 = rate_matrix.shape[1]
    # Initialization
    sparsity_num = 0
    sparsity_den = 0

    tall = np.sum(time_bin)

    # Occupancy probability
    pall = time_bin/tall

    for i in range(n1):
        for j in range(n2):

            # rate of each bin
            lx = rate_matrix[i, j]
            if probability == 'uniform':
                # Probability of occuping each bin - uniform
                p = 1.0/(n1*n2)
            elif probability == 'random':
                # Probability of occuping each bin - non-uniform
                p = pall[i, j]

            # numerator
            sparsity_num += p*lx
            # denominator
            sparsity_den += p*(lx**2)

    if sparsity_den == 0:
        sparsity_den = 1e-15

    sparsity = (sparsity_num)**2/float(sparsity_den)

    return sparsity


def sparsity_index2(rate_matrix):
    # Size of field
    n1 = rate_matrix.shape[0]
    n2 = rate_matrix.shape[1]
    # Initialization
    sparsity_num = 0
    sparsity_den = 0
    for i in range(n1):
        for j in range(n2):

            # rate of each bin
            lx = rate_matrix[i, j]
            # numerator
            sparsity_num += lx
            # denominator
            sparsity_den += lx**2

    if sparsity_den == 0:
        sparsity_den = 1e-15
    L = n1

    sparsity = 1 - (1.0/L)*(sparsity_num**2 /
                            float(sparsity_den))*(float(L)/(L-1))

    return sparsity


def peak_frequency(rate_matrix):
    return np.max(rate_matrix)


def field_size(rate_matrix, relfreq=1.0, track_length=200):
    if rate_matrix.shape[0] == 1 or rate_matrix.shape[1] == 1:
        peak = np.argmax(rate_matrix)

        # left from peak + peak
        counter1 = 0
        while rate_matrix[peak-counter1] >= relfreq:
            counter1 += 1
            if peak-counter1 < 0:
                counter1 -= 1
                break
        # right from peak
        counter2 = 0
        while rate_matrix[peak+counter2] >= relfreq:
            counter2 += 1
            if peak+counter2 >= rate_matrix.shape[0]:
                counter2 -= 1
                break

        size = counter1+counter2
        size -= 1

    else:
        pass

    if peak != 0 or peak != track_length:
        mean_in_place = np.mean(rate_matrix[peak-(counter1-1):peak+counter2+1])
    elif peak == 0:
        mean_in_place = np.mean(rate_matrix[peak:peak+counter2+1])
    elif peak == track_length:
        mean_in_place = np.mean(rate_matrix[peak-(counter1-1):])

    if peak - counter1 == 0:
        mean_out_place1 = 0
    else:
        mean_out_place1 = np.mean(rate_matrix[:peak-(counter1-1)])
    if peak + counter2+1 == track_length:
        mean_out_place2 = 0
    else:
        mean_out_place2 = np.mean(rate_matrix[peak+counter2+1:])

    mean_out_place = np.mean([mean_out_place1, mean_out_place2])
    return size, mean_in_place, mean_out_place


def upper_tri_indexing(A):
    '''
    A: input matrix
    returns its upper triangular elements
    excluding the diagonal.
    A should be square matrix
    return only upper tiangular elements, without diagonal
    '''

    if len(A.shape) != 2:
        sys.exit("Input is not two-dimensional.")

    if (A.shape[0] != A.shape[1]):
        sys.exit("Matrix is not square.")

    m = A.shape[0]
    r, c = np.triu_indices(m, 1)
    return A[r, c]


def stability_index(x, y=None):
    if y is None:
        A = upper_tri_indexing(np.corrcoef(x.squeeze(), rowvar=1))
    else:
        A = np.corrcoef(x.squeeze(), y.squeeze())[0][1]

    return A


def meanfilt(x, k):
    '''
    https://gist.github.com/bhawkins/3535131
    '''
    """Apply a length-k median filter to a 1D array x.
    Boundaries are extended by repeating endpoints.
    """
    assert k % 2 == 1, "Median filter length must be odd."
    assert x.ndim == 1, "Input must be one-dimensional."
    k2 = (k - 1) // 2  # left and right parts of filter
    # Padding array
    y = np.zeros((len(x), k), dtype=x.dtype)
    # Original array in the middle
    y[:, k2] = x
    for i in range(k2):
        j = k2 - i
        y[j:, i] = x[:-j]
        y[:j, i] = x[0]
        y[:-j, -(i+1)] = x[j:]
        y[-j:, -(i+1)] = x[-1]

    return np.mean(y, axis=1)


def spatial_coherence(A, window):
    '''
    spatial coherence adapted from 
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4068307/pdf/fnbeh-08-00222.pdf
    rate_matrix: rate map per cell
    window: number of neighbors +1 position
    '''
    A_sm = meanfilt(A, window)
    R = np.corrcoef(A, A_sm)[0][1]
#    Z = np.math.atanh(R)
    return R

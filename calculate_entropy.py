from numpy import exp, isnan, log
from numpy.linalg import eigh
from scipy.special import xlogy

def m_entropy(G):
    
    '''
        G -> the input graph.
    '''
    
    M = G.M
    eigenvalues = G.eigenvalues
    
    f = eigenvalues / (2*M)
    
    H = -xlogy(f, f)/log(2)
    H[isnan(H)] = 0
    
    return H

def b_entropy(G, beta = 1):
    
    '''
        G -> the input graph.
        beta -> a time-scale parameter
    '''
    
    eigenvalues = G.eigenvalues
    
    f = exp(-beta*eigenvalues)
    f /= f.sum()
    
    H = -xlogy(f, f)/log(2)
    H[isnan(H)] = 0
    
    return H
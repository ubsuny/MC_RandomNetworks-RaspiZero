from numpy import isnan
from numpy.linalg import eigh
from scipy.special import xlogy

def m_entropy(G):
    
    '''
        G -> the input graph.
    '''
    
    M = G.M
    eigenvalues = G.eigenvalues
    
    f = eigenvalues / (2*M)
    
    H = -xlogy(f, f)
    H[isnan(H)] = 0
    
    return H
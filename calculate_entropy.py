from numpy import argsort, argwhere, exp, isnan, log, zeros
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

def edge_rankings(G, beta = 1):
    
    Hs = zeros((G.edges.shape[0]))
    
    for idx, edge in enumerate(G.edges):
        
        jdx = argwhere(edge == G.edges)
        
        A = G.edge_removal(edge, jdx)
        
        G.update_laplacian()
        Hs[idx] = b_entropy(G, beta = beta).sum()
        G.edge_addition(edge, jdx, A)
        
    edge_sort = argsort(argsort(-Hs))
    
    return Hs, edge_sort
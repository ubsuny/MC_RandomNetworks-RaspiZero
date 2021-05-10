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

    '''
        G -> the graph object
        beta -> the time-scale parameter
    '''
    
    Hs = zeros((G.edges.shape[0]))
    
    for idx, edge in enumerate(G.edges.copy()):
        
        A = G.edge_removal(edge, 0)
        
        G.update_laplacian()
        Hs[idx] = b_entropy(G, beta = beta).sum()
        G.edge_addition(edge, Hs.shape[0] - 1, A)
        G.update_laplacian()
        
    edge_sort = argsort(-Hs)
    
    edge_ranks = argsort(edge_sort)
    
    return Hs, edge_sort, edge_ranks


def top_ranked(sorts, locs, cutoff = 0.10):

    '''
        sorts -> the way that the edges are sorted by entropy value
        locs -> the locations of the edges.  1 is community, 0 if connecting
        cutoff -> the percentile to cutoff the edge rankings
    '''
    
    cutoff_idx = int(cutoff*sorts.shape[0])
    if cutoff_idx == sorts.shape[0]:
        is_connecting = ~locs[sorts]
    else:
        is_connecting = ~locs[sorts][:cutoff_idx]
    
    is_connecting = is_connecting.sum()/is_connecting.shape[0]
    
    return is_connecting*100
    

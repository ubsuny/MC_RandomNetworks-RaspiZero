from numpy import append, argwhere, arange, array, diag, delete, ones, random, sum, triu
from matplotlib import pyplot as plt
import networkx as nx

class Erdos_Renyi_GNP:
    def __init__(self, N, p, A = None, self_edges = False):
        
        
        '''
            N -> the number of nodes, N >= 2, or vertices.  This must be declared on start-up to create a graph.
            p -> the probability on [0, 1] an edge forms between two nodes.  This must be declared on start-up.
            self_edges -> a Boolean whether or not a node can connect to itself.  In most applications, no.
        
        '''
        
        self.N = N
        self.p = p
        
        '''
            A -> the adjacency matrix, or the representation of all edge connections.  Index 0 refers to Node 0.
            D -> the degree matrix, or the sum of A of only one axes.  If 0, unconnected.  If N - 1, fully connected.
            L -> the Laplacian matrix as defined by L = D - A.  This allows for spectral analysis.
            M -> the total number of edges by summing A and then dividing by 2, since it mirrors across the diagonal.
        
        '''
        
        # This basically creates the adjacency matrix by creating an N x N random array
        # which has values from [0, 1].  So, it should be probabilistic enough to use < p
        # to remove all values greater than that.  Taking the upper triangular (triu)
        # allows us to force symmetry.  Although, this is one source of error
        # because it interrupts the probabilistic process.
        if A == None:
            self.A = triu(array(random.rand(N, N) < p, dtype = int))
            if self_edges == False:

                # If there are no self edges, we remove the diagonal
                # since entry [0, 0] for example refers to Node0 <-> Node0.
                self.A = self.A + self.A.T - 2*diag(diag(self.A))

            else:
                self.A = self.A + self.A.T - diag(diag(self.A))
        else:
            self.A = A
        
        # This sums the elements vertically then puts them in the diagonals
        # of a zero array.
        self.D = diag(sum(self.A, axis = 1))
        self.L = self.D - self.A
        
        # Since we assume this matrix is non-directional, then the lower and upper parts
        # have the same values, so we must divide by two after summing.
        self.M = sum(self.A)/2
        
        self.edges = argwhere(triu(self.A) != 0)
        self.potential_edges = argwhere((triu(1 - self.A) - diag(ones(self.N)) != 0))
     
    def plot_graph(self, figsize = (4, 4)):
        
        fig, ax = plt.subplots(1, 1, figsize = figsize)
        
        ax.spy(self.A)
        
        ax.set_ylabel('node ID, $y$', fontsize = 12)
        ax.set_xlabel('node ID, $x$', fontsize = 12)
        ax.xaxis.set_label_coords(0.5, 1.175)
        
        return fig, ax
    
    def plot_networkx(self, figsize = (4, 4), node_color = 'navy', node_alpha = 0.75, node_size = 100, edge_color = 'black', edge_alpha = 0.15):
        
        fig, ax = plt.subplots(1, 1, figsize = figsize)
        
        G = nx.Graph()
        
        G.add_nodes_from(arange(self.N))
        G.add_edges_from(self.edges)
        
        pos = nx.spring_layout(G)
        
        nx.draw_networkx_nodes(G, pos = pos, node_color = node_color, node_size = node_size, alpha = node_alpha, ax = ax)
        nx.draw_networkx_edges(G, pos = pos, edge_color = edge_color, alpha = edge_alpha, ax = ax)
        
        plt.axis('off')
        
        fig.tight_layout()
        
        return fig, ax
    
    def rewire_graph(self):
        
        # Choose at random since ER
        idx = random.randint(self.edges.shape[0])
        
        # Remove
        r_ij, r_ji = self.edges[idx]
        A_ij, A_ji = self.A[r_ij, r_ji].copy(), self.A[r_ji, r_ij].copy()
        
        self.A[r_ij, r_ji] = 0
        self.A[r_ji, r_ij] = 0
        
        self.potential_edges = append(self.potential_edges, self.edges[idx].reshape(1, 2), axis = 0)
        self.edges = delete(self.edges, idx, axis = 0)
        
        potential_idx = random.randint(self.potential_edges.shape[0])
        
        a_ij, a_ji = self.potential_edges[potential_idx]
 
        self.A[a_ij, a_ji] = A_ij
        self.A[a_ji, a_ij] = A_ji
        
        self.edges = append(self.edges, self.potential_edges[potential_idx].reshape(1, 2), axis = 0)
        self.potential_edges = delete(self.potential_edges, potential_idx, axis = 0)
        
        self.M = self.edges.shape[0]
    
from numpy import append, argwhere, arange, array, diag, delete, ones, random, sum, triu, zeros
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

        if A is None:
            self.A = triu(array(random.rand(N, N) < p, dtype = int))
            if self_edges == False:

                self.A = self.A + self.A.T - 2*diag(diag(self.A))

            else:
                self.A = self.A + self.A.T - diag(diag(self.A))
                
        else:  self.A = A

        self.D = diag(sum(self.A, axis = 1))
        self.L = self.D - self.A
        
        self.M = sum(self.A)/2
        
        self.edges = argwhere(triu(self.A) != 0)
        self.potential_edges = argwhere((triu(1 - self.A) - diag(ones(self.N)) != 0))
     
    def plot_graph(self, figsize = (4, 4)):

        '''

            figsize -> the dimensions in (x, y) of the graph.
            
        '''
        
        fig, ax = plt.subplots(1, 1, figsize = figsize)
        
        ax.spy(self.A)
        
        ax.set_ylabel('node ID, $y$', fontsize = 12)
        ax.set_xlabel('node ID, $x$', fontsize = 12)
        ax.xaxis.set_label_coords(0.5, 1.175)
        
        return fig, ax
    
    def plot_networkx(self, figsize = (4, 4), node_color = 'navy', node_alpha = 0.75, node_size = 100, edge_color = 'black', edge_alpha = 0.15):

        '''

            figsize -> the dimensions in (x, y) of the graph.
            node_color -> the color of the vertices of the networkx plot.
            node_alpha -> the transparency of the vertices, on a scale of [0, 1].
            node_size -> the size of the nodes.  This scales with figsize.
            edge_color -> the color the connecting edges.
            edge_alpha -> the transparency of the connecting edges, on a scale of [0, 1].
            
        '''
        
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
    
    def edge_removal(self, edge, idx):
        
        '''
            edge -> the edge in (i, j) to be removed.
            idx -> the location of that edge.         
        '''
        
        r_ij, r_ji = edge
        A_ij, A_ji = self.A[r_ij, r_ji].copy(), self.A[r_ji, r_ij].copy()
        
        self.A[r_ij, r_ji] = 0
        self.A[r_ji, r_ij] = 0
        
        self.potential_edges = append(self.potential_edges, edge.reshape(1, 2), axis = 0)
        self.edges = delete(self.edges, idx, axis = 0)
        
        return A_ij, A_ji
    
    def edge_addition(self, edge, idx, A):
        
        '''
            edge -> the edge in (i, j) to be added.
            idx -> the location of that edge.
            A -> the weights of the edge in (A_ij, A_ji) to be added.
        '''
        
        a_ij, a_ji = edge
 
        self.A[a_ij, a_ji] = A[0]
        self.A[a_ji, a_ij] = A[1]
        
        self.edges = append(self.edges, edge.reshape(1, 2), axis = 0)
        self.potential_edges = delete(self.potential_edges, idx, axis = 0)
        
    def rewire_graph(self):
        
        idx = random.randint(self.edges.shape[0])

        A_ij, A_ji = self.edge_removal(self.edges[idx], idx)
        
        potential_idx = random.randint(self.potential_edges.shape[0])
 
        self.edge_addition(self.potential_edges[potential_idx], potential_idx, (A_ij, A_ji))
    
        self.M = self.edges.shape[0]
    

    
class SBM(Erdos_Renyi_GNP):
    def __init__(self, N, p_in, p_out, k, A = None, self_edges = False, one_hot = None):
        
        '''
            p_in -> probability of edge forming inside communities.  Shared across all communities.
            p_out -> probability of edge forming within communities.  Shared across all communities.
            k -> number of communities.
            one_hot -> the one hot encoding for this particular community structure.
        '''
        
        self.p_in = p_in
        self.p_out = p_out
        self.N = N
        
        self.n = N//k
        if one_hot is None:
            self.one_hot = zeros((self.N, k))

            for dim in range(k):
                self.one_hot[(self.n*dim):(self.n*(dim + 1)), dim] = 1
                
        else: self.one_hot = one_hot
            
        self.in_mask = self.one_hot @ self.one_hot.T
        self.out_mask = 1 - self.in_mask
        
        self.A_in = Erdos_Renyi_GNP(self.N, self.p_in).A * self.in_mask
        self.A_out = Erdos_Renyi_GNP(self.N, self.p_out).A * self.out_mask
        
        self.inside_edges = argwhere(triu(self.A_in) != 0)
        self.outside_edges = argwhere(triu(self.A_out) != 0)
        
        self.potential_inside_edges = argwhere((triu(1 - self.A_in) * self.in_mask - diag(ones(self.N)) != 0))
        self.potential_outside_edges = argwhere((triu(1 - self.A_out) * self.out_mask != 0))
        
        self.A = self.A_in + self.A_out
        
        Erdos_Renyi_GNP.__init__(self, N, 0, A = self.A)
    
    def get_edge_colors(self, in_color, out_color):
        
        '''
            in_color -> the color of inside community edges
            out_color -> the color of between community edges.
        '''
        
        colors = array([])
        
        for edge in self.edges:
            i, j = edge
            if self.in_mask[i, j] == 1:
                colors = append(colors, in_color)
            else: colors = append(colors, out_color)
                
        return colors
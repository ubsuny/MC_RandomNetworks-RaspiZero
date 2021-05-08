import sys
sys.path.insert(1, '../')

from generate_network import *

# Initial conditions
N = 100
p = 0.25

# Generate graph
G = Erdos_Renyi_GNP(N, p)

def test_network_size():
    assert(G.A.shape == (N, N))

def test_edge_numbers():
    bar_m = lambda N, p: N*(N - 1)/2 * p
    std_m = lambda N, p: 0.10 * bar_m(N, p)

    assert((G.M < bar_m(N, p) + std_m(N, p)) and (G.M > bar_m(N, p) - std_m(N, p)))

    for _ in range(100):
        G1 = Erdos_Renyi_GNP(N, p)
        assert((G1.M < bar_m(N, p) + std_m(N, p)) and (G1.M > bar_m(N, p) - std_m(N, p)))

def test_total_edges():
    assert(G.edges.shape[0] == G.M)
    assert(G.edges.shape[0] + G.potential_edges.shape[0] == G.N*(G.N - 1)/2)

def test_rewiring():
    M = G.M
    for _ in range(100):
        G.rewire_graph()
        assert(G.M == M)
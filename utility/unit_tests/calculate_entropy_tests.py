import sys
sys.path.insert(1, '../')

from generate_network import *
from calculate_entropy import *
from numpy import log2

# Initial conditions
N = 100
p = 0.25

# Generate graph
G = Erdos_Renyi_GNP(N, p)

# M entropy tests:
def test_minimum_entropy():
    assert(m_entropy(G).sum() >= 0)

def test_maximum_entropy():
    assert(m_entropy(G).sum() <= log2(N))

def test_scaling_entropy():
    Ns = range(10, 100, 5)
    for idx, N in enumerate(Ns[1:]):
        G0 = Erdos_Renyi_GNP(Ns[idx], p)
        G1 = Erdos_Renyi_GNP(N, p)
        assert(m_entropy(G1).sum() >= m_entropy(G0).sum())

# B entropy tests:
def test_community_entropy():

    N = 100
    ks = range(1, 10, 2)
    p_in = 1
    p_out = 0

    for idx, k in enumerate(ks[1:]):
        G0 = SBM(N, p_in, p_out, k = ks[idx])
        G1 = SBM(N, p_in, p_out, k = k)

        assert(b_entropy(G0).sum() <= b_entropy(G1).sum())

def test_rewiring_entropy():

    N = 100
    k = 2
    p_in = 1
    p_out = 0

    G = SBM(N, p_in, p_out, k = k)
    H0 = b_entropy(G).sum()
    for t in range(1000):
        G.rewire_graph()
        if t % 25 == 0 and t != 0:
            H1 = b_entropy(G).sum()
            assert(H1 < H0)
            H0 = H1

test_rewiring_entropy()

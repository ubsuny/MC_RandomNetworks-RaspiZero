# Random Networks via Monte Carlo Simulations on the Raspberry Pi Zero W

### A Computational Physics II project by Jeremy Kazimer (jdkazime@buffalo.edu)

---

#### Introduction

Whether it's obvious or not, networks are everywhere.  Well, to be more precise, the information that they encode is everywhere.  What I mean by this is that networks on their highest level are composed of two properties:  nodes,  $N$, and edges, $M$.  More will be elaborated upon in the theory section, but essentially each node represents some object in a pre-defined space.  This could be a person, a neuron, or...well, anything, really.  Anything that takes on an identity and can be classified by said identity.   As for the other property, the edges, these encode the interactions between any set of nodes.  If, for example, two people follow each other on Twitter, then they have an edge between them.  These are the surface level properties of a network.  There are of course more, but these are the essentials.

Thus, taking a step back, it is more obvious that networks are everywhere.  To elaborate, with the right data as defined by information theory (CITE), a network can be constructed such that the nodes encode any set of objects and the edges the subsequent interactions.   To continue the Twitter analogy, take any $N$ twitter users.  These are the nodes, the objects, of our system.  Each of theses users, these nodes, have unique properties, such as their username, their profile picture, etc.  This is what it means for a node to encode information.  They don't have to encode information, but they inherently have the structure to do so.  Regardless, we can then define the existence of an edge on whether a user follows another user.  If they do and follow each other, then they share an edge.  If one person follows another, but not the other way around, then there exists only one edge from that user to the followed user.  These are directional edges.  There are also weighted edges such that perhaps an interaction between two nodes is more important than the other nodes.  Either way, if they do not follow each other, then there is no edge between that node.   

Really, what an edge represents here is the ability for information to spread; in a closed system, that is to say that there is no possibility of two completely separated nodes of communicating, information travels from node to node by their edges.  As a physical analogy, think of this like a set of interacting particles; if there is no collision, that is to say particles are not interacting with each other, then energy is not spread around the system, assuming the lack of physical entropy.   However, if they are interacting through collision, then this energy, the information, is spread via this transfer.  

Ultimately, with the above assessment, the goal of this project is to characterize these networks on a purely random level.  That is to say, for a network of size $N$ it is generated randomly with edge probability $p$, the probability of an edge between two nodes, so that both the configuration and the number of edges $M$ is random each time, within the probability constraints.   With that, we can create an analogy between this and that of actual particle interaction via information-theoretic entropic analysis, among other things.  In order to make definitive statements about these random networks, the goal of this assignment, I will employ the use of Monte Carlo methods so that we can validate known theory and also make observations of our own.   As a further constraint, this will also be created such that it can be reasonably ran on a Pi Zero W, which introduces its own margin of error.



#### Theory

In order to generate networks, we must first define all of its properties such that the properties are understood.  First, speaking of nodes we define

$\mathcal{V} = \{v_1, v_2, \dots , v_n\}$

for $n \in \mathbb{N}$ such that $\mathcal{V}$ is the set of all nodes and $v_1$ and $v_2$ represent the first and second nodes, respectively.  The number of nodes $N$ is then just the cardinality of $\mathcal{V}$ such that $N = |\mathcal{V}|$ .  Then, for a single-layer Erdős–Rényi graph $G_{NP}$ (CITE), that is to say a two-dimensional graph with no coupling constant $D_x$ (CITE), this is defined as the adjacency matrix $A$ such that 

$G_{NP} = A = \begin{bmatrix} A_{ij}\end{bmatrix}$

for $ 0 \leq i, j < N$ .  Then, the number of elements in the adjacency would then be

$|G_{NP}| = |A| = |N|\times |N|$

Note that this is not the number of edges, but rather the number of all entries for this adjacency matrix or, in computer science terms, the size of the array.  From this point onward, we'll only refer to $G_{NP}$ as $A$.   We can then define the values of $A_{ij}$ such that

$A_{ij} = \begin{cases} 1: (i, j) \in E \\ 0: (i, j) \notin E\end{cases}$

where $E$ is the set of all edges such that 

$E = \{(i, j): p_{ij} \leq p\}$

for some particular edge probability $p_{ij}$ and a total edge probability $p$.  Basically, each entry of $A$, $A_{ij}$, has an associated probability $p_{ij}$ such that if it's less than the fixed parameter probability $p$ then there exists an edge.  We can then define thte total number of edges as

$M = |E|$ 

Therefore, the number of potential edges $P$, that is to say areas where edges don't exist such that they could exist if under some stochastic process, can be represented as 

$P = \frac{|A|}{2} - M$

since $|A|$ counts edges twice.  
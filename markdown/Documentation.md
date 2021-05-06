# Random Networks via Monte Carlo Simulations on the Raspberry Pi Zero W

## A Computational Physics II project by Jeremy Kazimer (jdkazime@buffalo.edu)

### Introduction

Whether it's obvious or not, networks are everywhere.  Well, to be more precise, the information that they encode is everywhere.  What I mean by this is that networks on their highest level are composed of two properties:  nodes,  $N$, and edges, $M$.  More will be elaborated upon in the theory section, but essentially each node represents some object in a pre-defined space.  This could be a person, a neuron, or...well, anything, really [@physics].  Anything that takes on an identity and can be classified by said identity.   As for the other property, the edges, these encode the interactions between any set of nodes.  If, for example, two people follow each other on Twitter, then they have an edge between them.  These are the surface level properties of a network.  There are of course more, but these are the essentials.

Thus, taking a step back, it is more obvious that networks are everywhere.  To elaborate, with the right data as defined by information theory [@networkentropy], a network can be constructed such that the nodes encode any set of objects and the edges the subsequent interactions.   To continue the Twitter analogy, take any $N$ twitter users.  These are the nodes, the objects, of our system.  Each of theses users, these nodes, have unique properties, such as their username, their profile picture, etc.  This is what it means for a node to encode information.  They don't have to encode information, but they inherently have the structure to do so.  Regardless, we can then define the existence of an edge on whether a user follows another user.  If they do and follow each other, then they share an edge.  If one person follows another, but not the other way around, then there exists only one edge from that user to the followed user.  These are directional edges.  There are also weighted edges such that perhaps an interaction between two nodes is more important than the other nodes.  Either way, if they do not follow each other, then there is no edge between that node.   

Really, what an edge represents here is the ability for information to spread; in a closed system, that is to say that there is no possibility of two completely separated nodes of communicating, information travels from node to node by their edges.  As a physical analogy, think of this like a set of interacting particles [@dynamics]; if there is no collision, that is to say particles are not interacting with each other, then energy is not spread around the system, assuming the lack of physical entropy.   However, if they are interacting through collision, then this energy, the information, is spread via this transfer.  

Ultimately, with the above assessment, the goal of this project is to characterize these networks on a purely random level.  That is to say, for a network of size $N$ it is generated randomly with edge probability $p$, the probability of an edge between two nodes, so that both the configuration and the number of edges $M$ is random each time, within the probability constraints [@dane].   With that, we can create an analogy between this and that of actual particle interaction via information-theoretic entropic analysis, among other things.  In order to make definitive statements about these random networks, the goal of this assignment, I will employ the use of Monte Carlo methods so that we can validate known theory and also make observations of our own.   As a further constraint, this will also be created such that it can be reasonably ran on a Pi Zero W, which introduces its own margin of error.

### Limitations of the Raspberry Pi Zero W

The deployment of this application on the Raspberry Pi Zero W is not a particularly difficult feat; none of the involved libraries here are difficult to install or unable to be installed on the Raspberry Pi Zero W.  Nor are the calculations time intensive; most of these calculations are upper bounded by $\mathcal{O}(N^4)$ such that $N$ is the network size.

Really, the limiting factor here is space.  Since graph structures are represented in $y$ dimensional arrays such that $y \geq 2$, their space usage does not scale particularly well.  Notably, since most graphs are represented symmetrically, at least the ones we'll be dealing with, there will always be a spatial demand of at least $N \times N$.   If this is a multilayer structure [@multilayer], then this is brought into the third dimension by some number of layers $l$.  Further, since each float is at minimum 4 bytes and maximum 16 bytes, the space allocation $S$ scales $4 \cdot (N \times N \times l) \leq S \leq 16 \cdot (N \times N \times l)$.

As such, much of this experimentations for the Raspberry Pi Zero must be done with smaller networks, say 100 nodes.  Although this seems like a lot, physical systems are in the thousands to millions of nodes so we will have difficulty presenting any tangible system in our experimental environment.  

There are, however, experimental libraries in Python designed for the Raspberry Pi Zero W, such as VideoCore, to optimize libraries such as NumPy.  The issue is that it's rather experimental so any incorrect configuration could evaporate the board.  Not literally, but it would at worst brick the system.  Especially in such a complicated set of calculations, it's not worth the risk of trying to set it up.  With that, the solution we'll be utilizing is instead to just limit space usage.  And by virtue of working on smaller networks, the computation time will also speed up since it depends on network size.

With that, discussion on theory can begin.

### Network Theory

#### Graph Structure

In order to generate networks, we must first define all of its properties such that the properties are understood.  First, speaking of nodes we define

$\mathcal{V} = \{1, 2, \dots , n\}${#eq:description}

for $n \in \mathbb{N}$ such that $\mathcal{V}$ is the set of all nodes [@dane] and $1$ and $2$ represent the first and second nodes, for example.  The number of nodes $N$ is then just the cardinality of $\mathcal{V}$ such that $N = |\mathcal{V}|$ .  Then, for a single-layer graph Erdős–Rényi $G_{NP}$ , [@dane] that is to say a two-dimensional graph with no coupling constant $D_x$ [@multilayer], this is defined as the adjacency matrix $A$ such that 

$G_{NP}\left(\mathcal{V}, E\right) = A = \begin{bmatrix} A_{ij}\end{bmatrix}${#eq:description}

for $0 \leq i, j < N$ .  Then, the number of elements in the adjacency would then be

$|G_{NP}| = |A| = |N|\times |N|${#eq:description}

Note that this is not the number of edges, but rather the number of all entries for this adjacency matrix or, in computer science terms, the size of the array.  From this point onward, we'll only refer to $G_{NP}$ as $A$.   We can then define the values of $A_{ij}$ such that

$A_{ij} = \begin{cases} 1: (i, j) \in E \\ 0: (i, j) \notin E\end{cases}${#eq:description}

where $E$ is the set of all edges [@dane] such that 

$E = \{(i, j): p_{ij} \leq p\}${#eq:description}

for some particular edge probability $p_{ij}$ and a total edge probability $p$.  Basically, each entry of $A$, $A_{ij}$, has an associated probability $p_{ij}$ such that if it's less than the fixed parameter probability $p$ then there exists an edge.  We can then define the total number of edges as

$M = |E|${#eq:description}

For an unweighted and undirected matrix, this is the same as saying $2M = \sum_{ij} A_{ij}$, due to the symmetry of the adjacency matrix [@dane].  The expected number of edges then, for an $A$ graph, would be

$\bar{M} = \dfrac{N(N - 1)}{2} \cdot p${#eq:description}

because, without self edges, that is we disallow a node from being connected to itself, there are $N$ nodes.  However, we remove the locations where a node can connect itself such that there remains $N - 1$ spots remaining.  Then, we divide by $2$, due once again to the symmetry of an undirected and unweighted adjacency matrix [@avrim].  This leaves us with $\frac{N(N - 1)}{2}$ edges.  Since this is a purely probabilistic process, we multiply this term by $p$ since generally an expectation function takes on the form

$\bar{x} = x \cdot p${#eq:description}

We can then relate this to the number of potential edges $P$, that is to say areas where edges don't exist such that they could exist if under some stochastic process, can be represented as 

$P = \frac{N(N - 1)}{2} - M${#eq:description}

since $|A|$ counts edges twice.   We can then define the degree matrix, that is the number of edges connected to each node, as 

$D = \text{diag}[d_1, d_2, \dots d_N]${#eq:description}

where $d_i = \sum_i A_{ij}$ for the $i$-th node [@dane].  The unnormalized Laplacian matrix $L$ is then given by the equation

$L = D - A${#eq:description}

This is particularly useful to us in that this is the foundation of spectral analysis.  That is to say, we use the Laplacian in order to extract its eigenvalues $\lambda$ and eigenvectors $\vec{v}$ for a variety of uses.  This also appears in fields such as machine learning and computer vision, but for our uses this is simply the basis of spectral analysis.

#### Stochastic Block Model

A stochastic block model (SBM) is a structure of graphs that describe random graphs.   Here, it is structured in such a way that the diagonal blocks of the random graph are the communities, so to speak, of these graphs and the off diagonals represent the edges connecting these communities together.  (CITE)

More generally, we look at the case where all of the communities are equally sized; that is to say, the size of each community $n_k$ is defined by

$n_k = \frac{N}{k}${#eq:description}

where $k$ is the number of communities.  It must be true then that $\sum_i^k n_k = N$.  Although this describes the network size, this does not describe the number of edges.  Typically, we describe the probability structure of the graph by the matrix

$\begin{bmatrix} p_{in}^{(0)} & p_{out}^{(0, 1)} & \dots & p_{out}^{(k - 1, k)} \\ p_{out}^{(1, 0)} & p_{in}^{(1)} & \dots & \dots \\ \dots & \dots & \dots & \dots \\ p_{out}^{(k, k - 1)} & \dots & \dots & p_{in}^{(k - 1)} \end{bmatrix}${#eq:description}

Note that the sum of

$\sum_{i}^{k - 1} \sum_{j}^k p_{out}^{(i, j)} + \sum_i p_{in}^{(i)}$ {#eq:description}

does not have to sum to unity, $1$.  Here, $p_{out}$ describes the probability of edges forming between the communities or what we'll describe as the connecting edges.  Then, $p_{in}$ describes the edges within communities, or what we'll describe as the community edges.  We can then decide the total edges within a community and connecting a community separately such that

$\bar{M}_{in} = \sum_{i}^{k - 1} \dfrac{n_k (n_k + 1)}{2} \cdot p_{in}^{(k)}${#eq:description}

and

$\bar{M_{out}} = \dfrac{1}{2}\sum_i^{k - 1} \sum_{j}^k n_k\cdot p_{out}^{(i, j)}${#eq:description}

Here, each off-diagonal block is capable of having diagonal elements and our summation assumes that the graph is unweighted and undirected, so we can just divide by two to account for duplicate blocks (i.e., $p_{out}^{(0, 1)} = p_{out}^{(1, 0)}$, or more simply a symmetric system).  Of course, then, the total number of edges would then just be

$\bar{M} = \bar{M_{out}} + \bar{M_{in}}${#eq:description}

with the reminder that we look at the expectation here, since forming a $G_{NP}$ is inherently random.  Overall, an SBM can be viewed as the superposition of several different random graphs of equal network size, but varying edge probabilities.  More intuitively, an SBM is analogous to, say, a typical physical partition problem.  

Suppose each community represents a separate partition (CITE).  Each community edge then represents the potential that exists between two particles or nodes.  In a closed system, i.e. all $p_{out} = 0$, there is no energy exchange between the partitions, hence a lack of edges.  As such, as $p_{out}$ increases throughout overall the system, the more energy distribution that there is since the partitions are opened and now the separate communities are allowed to interact.

That is all that will be said for now.  Later, this will be elaborated upon.

#### Rewiring

AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

### Von Neumann Entropy

#### Spectral Analysis

To begin, spectral analysis is the examination of eigenvalues from some arbitrary system (CITE).   Since we are dealing with graph structure, any spectral analysis deals primarily with that of the eigenvalues of the Laplacian matrix $L$ (CITE).   There are two types of the Laplacian we typically handle, the unnormalized Laplacian $L$ and the normalized Laplacian $L^{sym}$.  

Here, we'll only speak of the unnormalized Laplacian.  There is nothing particularly wrong with the normalized, it's just that the unnormalized appears more commonly in nature (CITE).  Then, the eigenvalues $\{\lambda_i\}$ and eigenvectors $\{\vec{v}\}$ can be found arbitrarily with any typical algorithm.

In network theory, the eigenvalues by themselves are indicative of a few things.  The first, the edge probability and number of nodes.  This is because in a typical $G_{NP}$, the average eigenvalue floats around

$\bar{\lambda_i} \approx N\cdot p${#eq:description}

For an SBM, this is

$\bar{\lambda_i} = N \cdot p_{in}${#eq:description}

assuming that each $p_{in}^{(k)}$ is of equal value.  This is because $p_{out}$, assuming that $p_{out}^{(k - 1, k)}$ is the same for all $k$, is more indicative of the community structure than the edge probabilities and nodes.  

Speaking of which, we can determine, for a $G_{NP}$, how many communities there are by the number of zero eigenvalues (CITE).  Or more specifically

$k = |\{\lambda_i: \lambda_i = 0 \text{ for all } i \in N - 1\}|${#eq:description}

Besides that, we often study how eigenvalues change under rewiring.  If, for example, we were examining a two community SBM and randomly rewire it such that it becomes a one community SBM, we would anticipate that the number of zero eigenvalues decreases from $2$ to $1$.  (CITE)

If examined on a much smaller scale, eigenvalues just in general give insight to the changing structure of the network.  It must be noted that we don't typically examine eigenvectors; this deals more with principal component analysis (PCA) (CITE) than, say, network theory's branch of spectral analysis.  

Of course, these can be examined, but there is no established theory to do so in the direction that this project is headed.  And, admittedly, I don't know about them to develop theory on my own.  

Otherwise, spectral analysis tends to ascend its eigenvalues by going into a more scalable form - von Neumann Entropy.

#### Von Neumann Entropy

Von Neumann Entropy (VNE) is but one way to interpret eigenvalues (CITE).  It originally comes from the equation of Shannon entropy (CITE) such that

$H = -\sum_i p_i \log p_i${#eq:description}

where $H$ is our entropy and $p_i$ is the probability of some event happening.  From a quantum perspective, think of an example such that we measure the probability of the location of a particle.  If we know that it exists, then it must exist in one such location so that 

$\sum_i p_i = 1$ {#eq:description} 

However, VNE takes this one step further by defining the probabilities as the trace of some density matrix, Tr($\rho$) (CITE), such that $\rho$ is the density matrix.  The entropy is then

$H = - \text{Tr}(\rho)\log_2 \text{Tr}(\rho)${#eq:description}

where base 2 comes from the binary representation of information (CITE).  Recent publications have claimed that Tr$(\rho)$ is really a function of the eigenvalues (CITE) such that 

$\text{Tr}(\rho)_i = \dfrac{\lambda_i}{2M}${#eq:description}

which would create an entropic model such that

$H_M = -\sum_i \dfrac{\lambda_i}{2M}\log_2\dfrac{\lambda_i}{2M}${#eq:description}

We denote $H_M$ as such, since it is an entropic model that depends on the number of edges.  Note that it is true (CITE)

$\sum_i \dfrac{\lambda_i}{2M} = 1${#eq:description}

Since this comes from the trace of the density matrix, it is rather evident that this entropy serves its purposes as a means to convey the density of the Laplacian.  Ultimately, this is indicative of the edge probability.  

As such, in a totally disconnected system such that $M = 0$, the entropy is minimized at $H_M = 0$ because all eigenvalues would be zero such that we define $0\log_2 0 = 0$.  However, when $p = 1$, that is to say the maximum $M$, all of the eigenvalues take on the same value such that

$\forall_{i, j \in N - 1} \lambda_i = \lambda_j$ {#eq:description}

So, with any typical entropy model, when all of the probabilities are the same the entropy is maximized.  It is reasonable then to conclude that the usefulness of this model isn't particularly high; it really just indicates the density of any Laplacian, which can be observed regardless if this Laplacian is known.  As such, it is only useful when the Laplacian is not known, in the rare case where only the eigenvalues are known.  

Because of this, data scientists sought to find a better entropy model; one that has more application to physical systems.  Recently, the following function was validated as being correct [@manlio]:

$H_\beta = -\sum_i \dfrac{e^{-\beta \lambda_i}}{Z} \log_2 \dfrac{e^{-\beta \lambda_i}}{Z}$ {#eq:description}

Here, $Z$ is the partition function (CITE) from statistical mechanics such that

$Z = \sum_je^{-\beta \lambda_j}${#eq:description}

where $\beta$ is some time-scale parameter.  Because of this parameter, we denote $H_\beta$ as such.  Note that in statistical mechanics $\beta$ is defined as 

$\beta = \dfrac{1}{k_BT}${#eq:description}

where $k_B$ is the Boltzmann constant and $T$ is the system temperature.  However, here it refers to the time domain as opposed to the energy domain, especially since we're dealing with information theory here.   Via quantum mechanics, this makes sense since energy and time cannot be known simultaneously (CITE), so in a way this is the complement to the energy-based parameter.

Regardless, it is apparently that $\sum_i \dfrac{e^{-\beta \lambda_i}}{Z} = 1$, since $Z$ is but the total of all of the denominator.  This model of entropy in particular is great for community detection [@manlio].  That is to say, for a dense enough network with community structure, an SBM for example, the entropy approximates to the following form (CITE):

$H_\beta \approx \log_2 k$

For example, a dense $G_{NP}$ with no community structure will evaluate to zero entropy.  However, totally independent nodes, that is to say $k = N$, will evaluate to maximal entropy, or the maximum number of communities.  As such, it is apparent how this could be useful in community detection.  This is where the model breaks away from the older model; it can detect community effortlessly giving it purpose above determining matrix density.

Really, the question here remains about how $\beta$ effects this entropy.  Consider $\beta$ as such:  this parameter is indicative of how long the nodes get to interact through edges.  So, small $\beta$ means that there is less time for the nodes to interact and large $\beta$ means that there is plenty of time to interact.

Later, this will be useful in edge ranking, more specifically in determining the most important edges.   Nevertheless, this is the extent of the entropy model; it is rather new, so the technicalities of it have yet to be worked out.  Especially since this partition of network science isn't the most popular.

#### Edge Ranking

START THIS

For now, we can make one more statement on $\beta$:  when $\beta$ is infinitely small, our $H_\beta \approx H_M$.  We can use a Taylor expansion to show this:

$\beta \approx 0: Z = \sum_j e^{-\beta \lambda_j} \approx \sum_j -\beta \lambda_j = -\beta \sum_j \lambda_j = -\beta (2M)$ {#eq:description}

Note that the statement $\sum_j \lambda_j = 2M$ actually comes from our original model (CITE).  Further,

$\beta \approx 0: e^{-\beta \lambda_i} \approx -\beta \lambda_i$ {#eq:description}

Then,

$\beta \approx 0: H_\beta \approx -\sum_i\dfrac{-\beta\lambda_i}{Z}\log_2\dfrac{-\beta\lambda_i}{Z} \approx -\sum_i\dfrac{-\beta\lambda_i}{-\beta(2M)}\log_2\dfrac{-\beta\lambda_i}{-\beta(2M)} =-\sum_i \dfrac{\lambda_i}{2M}\log_2\dfrac{\lambda_i}{2M}${#eq:description}

Ultimately, this explains why small $\beta$ is unable to distinguish important edges; it no longer acts an entropic model for community detection, but rather the density of the network.  More technically, the eigenvalues that would normally have the most important value to the entropy model are now zeroed out, meaning that they all weigh the same.  The rich structure that we can observe from spectral analysis vanishes, basically.

### Using Monte Carlo as an Explorative Tool

...



### Potential Applications to Real Physics



### Results

### Comparison to other Methods

Although there aren't many other analytical methods that are similar to network theory, spectral theory [@manlio] in particular, the one that stands out the most is principal component analysis (PCA) [@pca].  This family of algorithms in particular examines the eigenvalues and eigenvectors so that a fuller picture of the data is formed.  Typically, this is used to reduce down data to a scalable form such that loss of information is minimized.

Network theory and, by extension, spectral theory is a bit better than PCA in some cases because it allows for a more complex relationship to be developed by the data.  That is to say, complex networks are able to encode very specific interactions between nodes, especially when multiple layers, a multilayer or multiplex network, are introduced.  Since PCA reduces these features down to only the eigenvectors and eigenvalues typically, many of these rich relationships are lost in translation.  

Spectral theory is particularly important, especially with the introduction of Shannon entropy via von Neumann Entropy (VNE)  [@dane].  This could be perhaps be treated as an extension of both PCA and network theory, since entropy takes on the form of its eigenvalues.  Regardless, VNE is far more common a tool of analysis in network theory than it is in PCA.  Really, it in itself is reducing data down to a single value, akin to PCA, but this is used rarely to more so describe the entire system, as opposed to its individual components.  Perturbation theory [@dane] is also useful in studying how the network changes under slight modifications to the network, something that PCA doesn't exactly cover.

Nevertheless, it is the amount of information that can be gleaned from network theory that sets it above PCA.  The complex systems can also be captured, but at the cost of runtime complexity.  It is true that in general PCA is faster than network theory, since its main bound is calculating eigenvalues.  Whereas, with network theory, there are many methods from modularity, to rewiring, to edge ranking, which all have significant runtime complexities.  However, their availability as an analysis tool is more important than runtime, in a larger scale.  Of course, this may not be the case for the Raspberry Pi Zero W, but initial conditions can be set such that runtime complexity isn't a major factor.  

### Conclusion

....

### References
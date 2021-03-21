# Random Networks via Monte Carlo Simulations on the Raspberry Pi Zero W

### A Computational Physics II project by Jeremy Kazimer (jdkazime@buffalo.edu)

---

#### Introduction

Whether it's obvious or not, networks are everywhere.  Well, to be more precise, the information that they encode is everywhere.  What I mean by this is that networks on their highest level are composed of two properties:  nodes,  $N$, and edges, $M$.  More will be elaborated upon in the theory section, but essentially each node represents some object in a pre-defined space.  This could be a person, a neuron, or...well, anything, really.  Anything that takes on an identity and can be classified by said identity.   As for the other property, the edges, these encode the interactions between any set of nodes.  If, for example, two people follow each other on Twitter, then they have an edge between them.  These are the surface level properties of a network.  There are of course more, but these are the essentials.

Thus, taking a step back, it is more obvious that networks are everywhere.  To elaborate, with the right data as defined by information theory (CITE), a network can be constructed such that the nodes encode any set of objects and the edges the subsequent interactions.   To continue the Twitter analogy, take any $N$ twitter users.  These are the nodes, the objects, of our system.  Each of theses users, these nodes, have unique properties, such as their username, their profile picture, etc.  This is what it means for a node to encode information.  They don't have to encode information, but they inherently have the structure to do so.  Regardless, we can then define the existence of an edge on whether a user follows another user.  If they do and follow each other, then they share an edge.  If one person follows another, but not the other way around, then there exists only one edge from that user to the followed user.  These are directional edges.  There are also weighted edges such that perhaps an interaction between two nodes is more important than the other nodes.  Either way, if they do not follow each other, then there is no edge between that node.   

Really, what an edge represents here is the ability for information to spread; in a closed system, that is to say that there is no possibility of two completely separated nodes of communicating, information travels from node to node by their edges.  As a physical analogy, think of this like a set of interacting particles; if there is no collision, that is to say particles are not interacting with each other, then energy is not spread around the system, assuming the lack of physical entropy.   However, if they are interacting through collision, then this energy, the information, is spread via this transfer.  

Ultimately, with the above assessment, the goal of this project is to characterize these networks on a purely random level.  That is to say, for a network of size $N$ it is generated randomly with edge probability $p$, the probability of an edge between two nodes, so that both the configuration and the number of edges $M$ is random each time, within the probability constraints.   With that, we can create an analogy between this and that of actual particle interaction via information-theoretic entropic analysis, among other things.  In order to make definitive statements about these random networks, the goal of this assignment, I will employ the use of Monte Carlo methods so that we can validate known theory and also make observations of our own.   As a further constraint, this will also be created such that it can be reasonably ran on a Pi Zero W, which introduces its own margin of error.



#### Theory


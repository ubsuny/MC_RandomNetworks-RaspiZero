# Describing Code for the Subsequent Projects

### Description

The purpose of this document is to explain all of the code found throughout the repository that isn't covered by the docstrings.  Namely, much of the `.py` functions lack explanations as they exist to streamline code throughout the different notebooks.  As such, each subsequent section here will represent a different `.py` file, including its unit testing.



#### generate_networks.py

The only item here, the class, is here to store all of the network properties and to include its basic methods, mostly that of rewiring.   Of course, the imports come first:

```python
from numpy import append, argwhere, arange, array, diag, delete, ones, random, sum, triu
from matplotlib import pyplot as plt
import networkx as nx
```

The purpose of importing `numpy` in portions like these is to reduce the time it takes to begin the notebooks.  Namely, importing the entirety of the library takes time and space, so reducing it down to the only essential functions works wonderfully.  `numpy` is used primarily for its vector operations and data structures.

Further, `matplotlib` has the same treatment, by only importing a subsection in `matplotlib.pyplot`.  There are many plotting functions, so trying to reduce this down to only a few would pose its own challenge.  The same can be said for `networkx`, so importing them in bulk takes time, but is a necessary sacrifice.  Of course, both of these are used for their plotting capabilities.

Now, examining the actual class:

```python
class Erdos_Renyi_GNP:
```

For our purposes, the graphs that we are using can be represented by an Erdos-Renyi graph, through some type of manipulation.  As such, this class will act as the super class to any extended classes.  So, it takes in no other classes since it is the parent.

##### init

Inside of this, once all of the initializations have been executed, there exists this portion

```python
if A == None:
            self.A = triu(array(random.rand(N, N) < p, dtype = int))
            if self_edges == False:

                self.A = self.A + self.A.T - 2*diag(diag(self.A))

            else:
                self.A = self.A + self.A.T - diag(diag(self.A))
                
else:  self.A = A
```

which sets the adjacency matrix.   Since the class allows for an adjacency matrix to be passed beforehand, namely useful if copying the methods, we check to see if this flag has been changed from its default `None` to some value.   If it is, we assign the set value to the adjacency matrix.  If not, we form our own.

Effectively, `random.rand(N, N)` creates a matrix of size $N \times N$, with values $A_{ij} \in [0, 1]$.   If $A_{ij} < p$, where $p$ is the edge probability, then that entry is replaced with a $1$.  Otherwise, it's set to $0$.  This is the result of `array(random.rand(N, N) < p, dtype = int)`.   We take the `triu`, or upper triangular, of this result so that we can recombine it on the basis of self edges.  Of course, when looking at something like the number of edges this introduces a form of error, but by doing this it forces symmetry.

Now, `self_edges` is a flag determining whether or not a node can be connected to itself.  In our case, it never will.  So, we add `self.A` with its transpose `self.A.T` and substract `2*diag(diag(self.A))`.   Note that `diag(self.A)` grabs the diagonal elements of `self.A` and `diag(diag(self.A))` places those in the diagonal of a matrix of zeros with size $N \times N$.   Since `triu` include the $1$'s along the diagonal, they're added twice and as such we subtract twice the diagonal.

Otherwise, we only subtract once due to the symmetry allowing for the diagonal to be a multiple of two, at least for an undirected and unweighted graph.  To get the actual diagonal matrix, we execute

```python
self.D = diag(sum(self.A, axis = 1))
self.L = self.D - self.A
```

Here, `sum(self.A, axis = 1))` sums the value of `self.A` across the second axis (since everything starts at 0).  This could be done on the first axis as well, but for convention we use the second.  This produces a vector of size $N$, which represents the diagonal values of `self.A`.  We then use `diag` to force this vector into the diagonal of a zero matrix, like previously.  Subtracting `self.D` and `self.A` is just the Laplacian `self.L`, as defined in the Documentation.

Last of the `init` function is declaring the edges and their locations:

```
self.M = sum(self.A)/2

self.edges = argwhere(triu(self.A) != 0)
self.potential_edges = argwhere((triu(1 - self.A) - diag(ones(self.N)) != 0))
```

Note that `self.M` is just its definition as defined in the Documentation.  The other variables, however, are a bit more complex.  Note that `self.edges` represents all of the places where there's an edge.  Namely, `triu(self.A) != 0` gives a Boolean array of all the locations in the upper triangular where the values are not zero.  This is because we define the lack of an edge as $0$, whereas the existence of an edge can be any value except for $0$.  So, `argwhere` gives us the location of these.

As for `self.potential_edges`, this applies only for unweighted and undirected edges.  Basically, if we invert `self.A`, that is to say `1 - self.A`, this will show us all the locations where edges can be placed.  We subtract a diagonal of ones again, since `triu(1 - self.A)` will have ones in its diagonal.  Then, we apply the same rhetoric as finding `self.edges`.  As such, all of the `init` variables have been explained.
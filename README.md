# Random Networks as described by Monte Carlo Simulations on the Raspberry Pi Zero W

### Description

The purpose of this repository is to primarily explore Erdős–Rényi graphs (random networks) through the use of Monte Carlo methods.  A detailed explanation of what precisely a random network is can be found in _Theory_ of the documentation, as well as the merit of a Monte Carlo simulation.

Otherwise, in order to run these simulations it is first necessary to create a random graph.  This is in the `generate_network.py` file, which contains the classes necessary to generate these objects.  A surface level dive into this `.py` file can be seen in the `network_analysis.ipynb` notebook.  Here, I show how to use the `.py` file, verify that the methods agree with numerical generalizations of a network, and then do some very facetious runtime complexity analysis.

For now, this is all that I have.  In the future, I plan to do a deeper dive into first rewiring, to establish it as another random process by which Monte Carlo is applicable.  After that, entropy must be introduced as a medium for randomness.  Then, once all the random processes are pinpointed, I would like to begin the actual Monte Carlo simulations so that we can ascertain network properties and compare them to known formulations, possibly deriving approximations.

In the end, I would like this to be entirely ran on the Pi Zero W.  Consider that networks take up both space and time, since they're N-dimensional structures, so minimizing the space-time complexity would be beneficial for other people, but also the current research I'm engaged in.  The work here also benefits my current research, as formalizing some functions allows for streamlined analysis.

#### Other Files

For now, there is only a `Code.md` and `Documentation.md`.  The former explains all of the code not inside the Jupyter Notebooks.  The latter explains the entire experiment, from the ground up.  In the future, there will be locations for figures and videos.

#### Installation

In order to use any of this code, all that is needed is an installation of Jupyter Notebook via the Anaconda distribution with a version of 3.5 or higher.  The primary libraries at play here are `numpy`, `matplotlib`, and `networkx`, so if these don't come with the installation of Jupyter Notebook then these should be installed at the newest version.

#### Running the Code

Although this is designed for the Raspberry Pi Zero W, this can be ran on other systems.  Once all of the necessary libraries are installed, all of the code can be run as is, in the `network_analysis.ipynb` notebook.

#### Video Link

https://www.youtube.com/watch?v=aFxdBRfqHU4

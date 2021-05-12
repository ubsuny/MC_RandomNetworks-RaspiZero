# Random Networks as described by Monte Carlo Simulations on the Raspberry Pi Zero W

### Description

The purpose of this repository is to primarily explore Erdős–Rényi graphs (random networks) and SBM graphs (structured random networks) through the use of Monte Carlo methods.  A detailed explanation of what precisely a random network is can be found in _Theory_ of the `documentation/Documentation.pdf`, as well as the merit of a Monte Carlo simulation.

Otherwise, in order to run these simulations it is first necessary to create a random graph.  This is in the `utility/generate_network.py` file, which contains the classes necessary to generate these objects.  A surface level dive into this `.py` file can be seen in the `network_analysis.ipynb` notebook.  Here, I show how to use the `.py` file, verify that the methods agree with numerical generalizations of a network, and then do some very facetious runtime complexity analysis.

Once the random graph has been created, we can run some spectral analysis via the `utility/calculate_entropy.py` file in the `spectral_analysis.ipynb`.  This gives us proper insight into the structure of the graphs on a more fundamental level, transcending the input parameters.  Mostly, here you'll find eigenvalue analysis.  

After entropy has been established, we can use this utility file once again in `edge_rankings.ipynb` - here, we use entropy as a notion to define edge ranking.  We quickly observe what it means to rank an edge by its entropy and the implications it has for SBM's, in particular.  

Accompanying all of this is the `monte_carlo_simulations.ipynb`.  This notebook is primarily used to validate that Monte Carlo can, in fact, be used with random networks.  And, of course, how we can utilize it.  It's a very facetious look into its applications, its main goal in establishing a channel to use it.

Lastly, the bread and butter of the experiment: `toy_model.ipynb`.  This establishes a hypothetical scenario in which we can combine the above four notebooks into one compact examination.  Although it doesn't go much in depth, it is a toy model after all.  My main purpose in doing this is showing that an experimental model actually does agree with our expectations.  With that, we can look to the future and start applying this on large-scale, real data, if we would want a more practical application.

In the end, I would like this to be entirely ran on the Pi Zero W.  Consider that networks take up both space and time, since they're N-dimensional structures, so minimizing the space-time complexity would be beneficial for other people, but also the current research I'm engaged in.  The work here also benefits my current research, as formalizing some functions allows for streamlined analysis.

#### Other Files

As for other files, there is a `markdown/Code.md` and `markdown/Documentation.md`.  The former explains all of the code not inside the Jupyter Notebooks.  The latter explains the entire experiment, from the ground up.  In this same `markdown` folder, there is also the files used to convert these from a `.md` to a `.pdf`.  The respective `.pdf`'s can be found in the `documentation` folder.

There are also figures that I found rather important to the experiment in `figures`.  They don't have a description, but their titles are indicative of what they represent.  I think it's useful to have them there, at least, so that the person doesn't necessarily have to run the code every time to view the results.

Lastly, we also have very scarce unit tests in `utility/unit_tests`.  Since this application is random in nature, it is difficult to focus in on specific unit tests, so I only have a few.  In the future, I would like to use Monte Carlo to generate more - this would really highlight the usefulness of the methods.

#### Installation

In order to use any of this code, all that is needed is an installation of Jupyter Notebook via the Anaconda distribution with a version of 3.5 or higher.  The primary libraries at play here are `numpy`, `matplotlib`, and `networkx`, so if these don't come with the installation of Jupyter Notebook then these should be installed at the newest version.

#### Running the Code

Although this is designed for the Raspberry Pi Zero W, this can be ran on other systems.  Once all of the necessary libraries are installed, all of the code can be run as is, in all of the notebooks.

#### Video Link

Midterm: https://www.youtube.com/watch?v=aFxdBRfqHU4

Final: for the final submission.
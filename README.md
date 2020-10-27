# NetDynamics

This program is based on [RBNLab](http://turing.iimas.unam.mx/~cgg/rbn/index.html)

It let to explore some properties of different types of *Random Boolean Networks* (RBNs), specialy to look for degeneracy effects in RBNs. 

The `config` file defines the properties of the random Boolean network as listed below:
  - *nm*: Number of nodes (N) per module
    - when `m = 1`, the total number of nodes in the network is N
    - when `m > 1`, the total number of nodes in the network is N * m
  - *m*: Number of modules
  - *k*: Connectivity (0, 1]
  - *l*: Module connectivity level (0, 1]
  - *p*: Probability of 1's within rules (0, 1)
  - *topology*: The topology for the connection in the network:
    - Homogeneous = 0
    - Normal = 1
    - Scale free = 2
    - Coupled = 3
  - *redNodes*: Number of redundant Nodes.
  - *fd*: Function degeneracy degree [0, 1]
  - *id*: Input degeneracy degree [0, 1]
  - *od*: Output degeneracy degree [0, 1]
  - *iterations*: Nuber of nets to sample
  - *file*: Output file name with the network state, network connection
  - *operations*: topology,get_attractors,get_phase
  - *maxStates*: Maximum number of states to inspect an attractor before to look for another attractor (could be enormous and need a huge quantity of memory storing attractors)


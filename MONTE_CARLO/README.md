# MONTE CARLO
This directory uses the core algorithm from **../ISING** (with a few small changes) in order to calculate a Transition Probability Matrix (TPM) and the Effective Information of a given model

## Download and Run
In order to use the code:

* **Download:** the primary files main.cpp, functions.cpp, EI.cpp, and header.h (or the entire directory)
* **Compile:** the code must be compiled to make an executable
  * g++ -std=c++11 -O3 main.cpp functions.cpp EI.cpp -o run_sim.exe
* **Run:** the code can be run from the command line without any arguments
  * ./run_sim.exe


## Output
* **params.txt:** contains simulation parameters such as lattice size and number of steps
* **states.txt:** contains the binary state of the spin lattice as a function of time 
* **compressed_states.txt:** if **compression_flag = true**, this contains the compressed representation of state space (i.e. only the observed states are assigned labels). Else, this will be identical to states.txt
* **mapping.txt:** contains the map between integer and compressed representation of the state space
* **time_series.txt:** spin lattice as a function of time. This file can get very large so set **output_spins = false** for large sims


#### DON'T FORGET TO LINK EI.CPP!


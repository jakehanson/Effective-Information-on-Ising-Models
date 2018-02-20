# SEMI-ANALYTIC
This directory contains codes that calculate the Transition Probability Matrix (TPM) and Effective Information (EI) of an Ising Model without actually simulating a Monte Carlo walk through state space. Instead, this code analytically calculates transition probabilities using the energy difference between states and the Metropolis Update Rule.

## Download and Run
In order to use the code:

* **Download:** the primary files main.cpp, functions.cpp, EI.cpp, and header.h (or the entire directory)
* **Compile:** the code must be compiled to make an executable
  * g++ -std=c++11 -O3 main.cpp functions.cpp EI.cpp -o run_sim.exe
* **Run:** the code can be run from the command line without any arguments
  * ./run_sim.exe


## Output
* **TPM_file.txt:** contains the TPM resulting from the semi-analytical calculations
* **EI_file.txt:** contains EI and model parameters

#### DON'T FORGET TO LINK EI.CPP!


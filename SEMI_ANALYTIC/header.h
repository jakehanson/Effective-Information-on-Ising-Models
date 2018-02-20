#pragma once
#include <vector>
#include <iostream>
#include <random>
#include <stdio.h>
#include <map>
#include <fstream>


std::tuple<int, int, int, int> get_nn(int index,int n_rows,int n_cols);

/* Function to calculate transition probability using metropolis hastings */
void update_TPM(long int N_states, double T, long int index1, long int index2, double diff, std::vector<std::vector<double>> &TPM);

/* Function to calculate Effective Information (Hoel, 2017) */
double get_EI(std::vector<double> ID,std::vector<double> ED,std::vector<std::vector<double>> TPM);

/* Function to calculate energy of a spin state */
template <size_t bitsetsize>
double get_E(std::bitset<bitsetsize> state){
	
	int n_rows = sqrt(state.size());
	int n_cols = sqrt(state.size());

	double E_tot; // energy value
	int left_nn,right_nn,up_nn,down_nn; // nearest neighbor values
	int left_index,right_index,up_index,down_index; // nearest neighbor indices

	/* For each state find the local energy and add it to the total */
	E_tot = 0;
	double nn_sum = 0;
	double value = 0; // holds the value of the spin matrix at the location of the index in question
	for(int index=0;index<state.size();index++){
		std::tie(left_index, right_index, up_index, down_index) = get_nn(index,n_rows,n_cols); // gets nn indices
		
		//switch 0's to -1's
		if(state[left_index]==0){
			left_nn = -1;
		}else{
			left_nn = 1;
		}
		if(state[right_index]==0){
			right_nn = -1;
		}else{
			right_nn = 1;
		}	
		if(state[up_index]==0){
			up_nn = -1;
		}else{
			up_nn = 1;
		}
		if(state[down_index]==0){
			down_nn = -1;
		}else{
			down_nn = 1;
		}
		if(state[index]==0){
			value = -1;
		}else{
			value = 1;
		}

		// calculate energy sum over nearest neighbors
		nn_sum = value*(left_nn+right_nn+up_nn+down_nn);
		E_tot = E_tot - 1/2.*nn_sum; // total energy of this state
	}

	return E_tot;
}

#include "header.h"
#include <cmath>
#include <stdexcept>

/* Function to tick transitions and return the Transition Probability Matrix */
std::vector<std::vector<double>> get_TPM(std::vector<unsigned long long> state_array,long long int size, long N_steps){

	std::vector<std::vector<double>> TPM;
	unsigned long long current_state;
	unsigned long long next_state;
	unsigned long long N_states;

	// Initialize
	TPM = std::vector<std::vector<double>>(size,std::vector<double>(size,0.0));
	N_states = state_array.size();

	// Tick Transitions, skipping transitions between sims
	for(int i=0;i<N_states-1;i++){
		if((i+1)%N_steps != 0 or i == 0){
			current_state = state_array[i];
			next_state = state_array[i+1];
			TPM[current_state][next_state]++;			
		}
	}

	return TPM;
}


/* Function to calculate the Effective Information */
double get_EI(std::vector<double> ID,std::vector<double> ED,std::vector<std::vector<double>> TPM){
	
	double sum = 0.;
	for(int i = 0;i<TPM.size();i++){
		double D_KL = 0.;
		for(int j = 0;j<TPM.size();j++){
			if(TPM[i][j] != 0.){
				D_KL += TPM[i][j]*std::log2(TPM[i][j]/ED[j]);
			}
		}
		sum = sum + D_KL*ID[i];
	}

	return sum;
}



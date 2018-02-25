#include "header.h"

// This code calculates pairwise energy between states sequentially
// It does not randomly draw the next state
// Monte Carlo Metropolis algorithm should converge to these TPMs
// We could cut time by roughly half by only calculating upper diagonal for TPM


/* Main */
int main(int argc, char** argv)
{
	int n_rows = 3; // number of rows
	int n_cols = 3; // number of columns
	double T_c = 2.26918531421; // critical temp of ising model (coupling = +1, dimless temperature)
	double T = 10*T_c;  // Temperature of the ising model
	const int N_spins = 9; // total number of spins (MUST BE INITIALIZED MANUALLY)

	try{
		if(N_spins != n_rows*n_cols){
			throw std::runtime_error("N_spins must be equal to n_rows*n_cols.");
		}
		if(n_rows != n_cols){
			throw std::runtime_error("Rows and columns must match.");
		}
		if(n_rows > 4 or n_cols > 4){
			throw std::runtime_error("State space is too large!");
		}
	}catch(std::runtime_error &e){
		std::cerr << "RUNTIME ERROR: " << e.what() << std::endl;
		return 1;
	}

	long int N_states = pow(2,n_rows*n_cols);
	std::vector<std::vector<double>> TPM(N_states,std::vector<double>(N_states,0.)); // 2d array for the transition probability matrix

	/* Open file for TPM */
	std::ofstream TPM_file("TPM_file.txt");
	TPM_file << "Rows = " << n_rows << std::endl;
	TPM_file << "Cols = " << n_cols << std::endl;
	TPM_file << "T = " << T << std::endl;
	TPM_file << "TPM:" << std::endl;

	/* Open file for EI and Params */
	std::ofstream EI_file("EI_file.txt");
	EI_file << "Rows" << "\t" << "Cols" << "\t" << "T" << "\t" << "EI" << std::endl;

	/* Calculate Delta E and TPM*/
	std::cout << "Starting.." << std::endl;
	for(long int i=0;i<N_states;i++){
		std::bitset<N_spins> state1(i); // initializes a bitset in the form of a binary string representing integer i
		// Flip each bit and get energy difference
		for(int j=0;j<N_spins;j++){
			std::bitset<N_spins> state2(i); // make a copy of state1
			state2.flip(j); // flip a bit
			// calculate energy difference
			double diff = get_E(state2)-get_E(state1); // energy difference
			update_TPM(N_states,N_spins,T,i,state2.to_ulong(),diff,TPM); // store in matrix
		}
		// Calculate probability of transition out of state
		double trans_prob = 0.;
		for(long int f=0;f<N_states;f++){
			trans_prob = trans_prob + TPM[i][f]; // sum transition probabilities
		}
		TPM[i][i] = 1.-trans_prob; // prob of staying in the same state
	}

	// Write TPM to file
	for(long int i=0;i<N_states;i++){
		double sum = 0.;
		for(long int j=0;j<N_states;j++){
			if(TPM[i][j] > 1 or TPM[i][j] < -1e-12){
				std::cout << "ERROR INVALID TRANSITION PROBABILITY: " << TPM[i][j] << std::endl;
			}else{
				if(j != N_states - 1){
					TPM_file << TPM[i][j] << "\t";
				}else{
					TPM_file << TPM[i][j] << "\n";
				}
				sum += TPM[i][j];
			}
		}
		if(sum - 1.0 > 1e-10){
			std::cout << "ERROR INVLALID TRANSITION PROBABILITY SUM: " << sum << std::endl;
		}

	}

	// Create Intervention Distribution (assuming ID=H_max)
	std::vector<double> ID(N_states, 0.0);
	for(int i=0;i<N_states;i++){
		ID[i] = 1./N_states; // append 1/n for each value of ID
	}

	// Get Effect Distribution ED=ID*TPM
	std::vector<double> ED(N_states, 0.0);
	for(int i=0;i<N_states;i++){
		double sum = 0;
		for(int j=0;j<N_states;j++){
			sum = sum + TPM[j][i]*ID[j];
		}
		ED[i] = sum;
	}

	// Get EI and write to file
	double EI;
	EI = get_EI(ID,ED,TPM);
	EI_file << n_rows << "\t" << n_cols << "\t" << T << "\t" << EI << std::endl;
	std::cout << "EFFECTIVE INFO = " << get_EI(ID,ED,TPM) << std::endl;

	// Close files
	TPM_file.close();
	EI_file.close();
	std::cout << "Done." << std::endl;
}



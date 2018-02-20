
#include "header.h"
#include <fstream>


/* Main Function to Evolve Ising Model and Return list of States */
int main(int argc, char** argv)
{

	// First get lattice size and temperature
	int n_rows;
	int n_cols;
	double T;
	double N_steps;
	long N_sims;
	long N_burns;

	if(argc == 1){
		n_rows = 3;
		n_cols = 3;
		T = 2.0;
		N_steps = 1e7;
		N_sims = 1;
		N_burns = 0;
	}
	if(argc == 6){
		n_rows = atoi(argv[1]);
		n_cols = atoi(argv[1]);
		T = atof(argv[2]);
		N_steps = atof(argv[3]);
		N_burns = atoi(argv[4]);
		N_sims = atoi(argv[5]);		
	}

	std::map<unsigned long long,unsigned long long> compression_map; // map storing minimal representation of state space
	std::vector<unsigned long long> state_array; // holds the state trajectory

	/* Make Sure Input Params are valid */
	try{
		if(argc != 6 and argc != 1){
			throw std::runtime_error("INVALID SYNTAX!");
		}		
		if(n_rows != n_cols){
			throw std::runtime_error("SQUARE LATTICES ONLY - N_ROWS MUST MATCH N_COLS. CHECK VALUES IN MAIN.CPP");  // raise error if lattice isn't square
		}
		if(T <= 0.){
			throw std::runtime_error("TEMPERATURE MUST BE POSITIVE!");
		}
		if(N_burns < 0){
			throw std::runtime_error("NUMBER OF BURNS MUST BE NON-NEGATIVE");
		}
		if(N_steps < N_burns){
			throw std::runtime_error("NUMBER OF STEPS MUST BE GREATER THAN NUMBER OF BURNS");
		}
	}catch(std::runtime_error &e){
		std::cerr << "RUNTIME ERROR: " << e.what() << std::endl;
		return 1;
	}


	for(int n=0;n<N_sims;n++){

		// Initialize
		double U4 = 0;
		double M4_avg = 0;
		double M2_avg = 0;
		double M_avg = 0;
		Ising_Model model = Ising_Model(n_rows,n_cols,T,N_steps,N_burns);

		// Evolve
		model.evolve();
		M4_avg = model.M4_tot/(model.current_step-model.num_burns);
		M2_avg = model.M2_tot/(model.current_step-model.num_burns);
		M_avg = model.M_tot/(model.current_step-model.num_burns);
		U4 = 1. - M4_avg/(3.*M2_avg*M2_avg);

		// Write results
		std::cout << n_cols << "\t" << T << "\t" << n << "\t" << M_avg << "\t" << U4 << std::endl;
	
	}

}

#include "header.h"
#include <fstream>

// LOOK INTO EXP CALL. COMPARE EXP TO UNIQUE ENERGY DIFFERENCES.

/* Main Function to Evolve Ising Model and Return list of States */
int main(int argc, char** argv)
{

	int n_rows = 8; // number of rows
	int n_cols = 8; // number of columns
	double T_c = 2.26918531421; // critical temp of ising model (coupling = +1, dimless temperature)
	double T = 0.1*T_c;  // Temperature of the ising model
	long N_steps = 1e6;  // Number of steps to take in a given simulation
	long N_sims = 10; // Number of times to reinitialize and run a simulation
	bool output_spins = false; // flag to output spin matrix at each time step

	/* If user specified input, use that instead */
	if(argc == 3){
		n_rows = atoi(argv[1]);
		n_cols = atoi(argv[1]);
		T = atof(argv[2]);
	}

	std::map<unsigned long long,unsigned long long> compression_map; // map storing minimal representation of state space
	std::vector<unsigned long long> state_array; // holds the state trajectory

	/* Make Sure Input Params are valid */
	try{
		if(n_rows != n_cols){
			throw std::runtime_error("SQUARE LATTICES ONLY - N_ROWS MUST MATCH N_COLS. CHECK VALUES IN MAIN.CPP");  // raise error if lattice isn't square
		}
		if(n_rows > 8){
			throw std::runtime_error("N_ROWS MUST BE LESS THAN 8 TO STORE STATES!");  // raise error if we have a lattice bigger than 8x8
		}
		if(T <= 0.){
			throw std::runtime_error("TEMPERATURE MUST BE POSITIVE!");
		}
		if(output_spins){
			double series_size = n_rows*n_cols*sizeof(int)*N_sims*N_steps/8/1e9; // memory storage required to store time series
			if(series_size > 1){
				std::cout << "WARNING!! EXTREMELY LARGE TIME SERIES FILE (" << series_size << " GB)\n\t CONSIDER SETTING OUTPUT_SPINS = FALSE" << std::endl;
			}
		}
		if(argc != 1 and argc != 3){
			throw std::runtime_error("INVALID ARGUMENTS! SYNTAX SHOULD BE: ./runsim.exe n T");
		}
	}catch(std::runtime_error &e){
		std::cerr << "RUNTIME ERROR: " << e.what() << std::endl;
		return 1;
	}

	std::cout << "Starting...\n";

	/* Store parameters in file */
	std::ofstream params("params.txt");
	params << "n_rows\tn_cols\tTemperature\tSteps\t\n";
	params << n_rows << "\t" << n_cols << "\t" << T << "\t" << N_steps<< "\n";
	params.close();


	/* Initialize and run simulations */
	//std::ofstream states("states.txt");
	//std::ofstream compressed_states("compressed_states.txt");
	std::ofstream time_series("time_series.txt");

	std::cout << "N_rows:\t" << n_rows << std::endl;
	std::cout << "N_cols:\t" << n_cols << std::endl;
	std::cout << "T:\t" << T << std::endl;
	std::cout << "Size of state space: " << pow(2.,n_rows*n_cols) << std::endl;\
	std::cout << "Steps\tStates\n"; // write some general info to screen

	state_array.reserve(N_steps*N_sims);
	for(int n=0;n<N_sims;n++){

		Ising_Model model = Ising_Model(n_rows,n_cols,T,N_steps); // Initialize new model at T

		/* Run Sim and Store Spin Matrices in File */
		if(output_spins == true){
			time_series << model.spin_matrix; // save first spin state to file
			model.evolve(time_series); // evolve it through N_steps
		}else{
			model.evolve();
		}

		// Update our compression map	
		compression_map = update_map(model.states,compression_map);

		// Write states of current model to file and update state array
		for(long j=0;j<N_steps;j++){
			//states << model.states[j] << std::endl;
			//compressed_states << compression_map[model.states[j]] << std::endl;
			state_array.push_back(compression_map[model.states[j]]); // add the compressed state to the state array
		}

		std::cout << (n+1)*N_steps << "\t" << compression_map.size() << std::endl;
		std::cout << "State Array Size\t" << state_array.size() << std::endl;

	
	}

	/* Close open files */
	//states.close();
	//compressed_states.close();
	time_series.close();

	/* Write compression map to file */
	// std::ofstream mapping("mapping.txt");
	// for(std::map<unsigned long long,unsigned long long>::iterator it=compression_map.begin(); it!=compression_map.end(); ++it){
 //    	mapping << it->first << " => " << it->second << '\n';
	// }
	// mapping.close();
	std::cout << "Done.\n"; // End Main

}
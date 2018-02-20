#include "header.h"
#include <cmath>
#include <stdexcept>

// This version of Metropolis update on Ising has a more intuitive notation than original for Delta E and nearest neighbor states

/* Define our constructor. Pass Args to this to create an instance of Ising Model Class */
/* Note: Its a constructor because there's no return type and method name 'Ising_Model' matches class name. */
Ising_Model::Ising_Model(int n_rows, int n_cols, double T, long N_steps,long N_burns){

	current_step = 0; // index for current step
	spin_matrix = std::vector<std::vector<int>>(n_rows,std::vector<int>(n_cols,0)); // initialize 2d matrix with null spins
	num_rows = n_rows;
	num_cols = n_cols;
	num_steps = N_steps;
	num_burns = N_burns;
	temp = T;
	M4_tot = 0.;
	M2_tot = 0.;
	M_tot = 0.;

	/* Initialize Spin Matrix with +1 and -1 randomly drawn */
	std::random_device rd; // used to obtain seed for random number engine
	std::mt19937 gen(rd()); // standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> uni_dis(0,1);
	for(int i=0;i<num_rows;i++){
		for(int j=0;j<num_cols;j++){
			int rand_int = uni_dis(gen);         // draw random int either 0 or 1
			if(rand_int == 1){
				 spin_matrix[i][j]= 1;
			}else{
				spin_matrix[i][j] = -1;  // if rand_int was 0 change it to -1
			}
		}
	}

	// If we don't have a burn in, calculate initialization
	if(num_burns == 0){
		double M = 0;
		for(int i=0;i<num_rows;i++){
			for(int j=0;j<num_rows;j++){
				M += spin_matrix[i][j];
			}
		}
		M = std::abs(M); // we want absolute value of magnetization
		M4_tot += pow(M,4.);
		M2_tot += pow(M,2.);
		M_tot += M/(num_rows*num_rows);
	}

}


// This method evolves and writes spins to file
void Ising_Model::evolve(void){
	int delta_E = 0;  // integer to hold the local energy value
	int row_index,col_index;
	int left_nn,right_nn,down_nn,up_nn; // integers to hold the values of the nearest neighbors

	std::random_device rd; // used to obtain seed for random number engine
	std::mt19937 gen(rd()); // standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> real_dis(0,1); // real distribution to draw from (0 to 1)
	std::uniform_int_distribution<> uni_dis(0,num_rows-1); // uniform distribution to draw from (0 to n_rows-1 inclusive)

	//For each step get the spin matrix and convert binary string to integer representation 
	for(int i=1;i<num_steps;i++){
		//time_series << i << "\n";
		current_step = i;

		/* Choose a node at random */
		col_index = uni_dis(gen); // random int between 0 and n_cols-1
		row_index = uni_dis(gen); // random int between 0 and n_rows-1

		// std::cout << "Run = " << i << std::endl;
		// std::cout << "\t Row = " << row_index << std::endl;
		// std::cout << "\t Col = " << col_index << std::endl;

		/* Get value of spin matrix at nearest neighbor sites using PBC*/
		if((row_index-1) >= 0){
			up_nn = spin_matrix[row_index-1][col_index];
			//std::cout << "\t Up = " << row_index-1 << " " << col_index << std::endl;
		}else{
			up_nn = spin_matrix[num_rows-1][col_index];
			//std::cout << "\t Up = " << num_rows-1 << " " << col_index << std::endl;
		}
		
		down_nn = spin_matrix[(row_index+1)%num_rows][col_index];
		//std::cout << "\t Down = " << (row_index+1)%num_rows << " " << col_index << std::endl;

		right_nn = spin_matrix[row_index][(col_index+1)%num_cols];
		//std::cout << "\t Right = " << row_index << " " << (col_index+1)%num_cols << std::endl;
		
		if((col_index-1) >= 0){
			left_nn = spin_matrix[row_index][col_index-1];
			//std::cout << "\t Left = " << row_index << " " << col_index-1 << std::endl;
		}else{
			left_nn = spin_matrix[row_index][num_cols-1];
			//std::cout << "\t Left = " << row_index << " " << num_cols-1 << std::endl;
		}

		/* Calculate energy and decide whether or not to flip (Metropolis Algorithm) */
		delta_E = 2*spin_matrix[row_index][col_index]*(left_nn+right_nn+down_nn+up_nn);
		if(delta_E < 0){
			spin_matrix[row_index][col_index] = -1*spin_matrix[row_index][col_index];
		}else{
			double r = real_dis(gen);
			if(r <= std::exp(-double(delta_E)/temp)){
				spin_matrix[row_index][col_index] = -1*spin_matrix[row_index][col_index];
			}
		}

		// Sum up spins to get magnetization of state
		if(i>num_burns){
			double M = 0;
			for(int i=0;i<num_rows;i++){
				for(int j=0;j<num_rows;j++){
					M += spin_matrix[i][j];
				}
			}
			M = std::abs(M); // we want absolute value of magnetization
			M4_tot += pow(M,4.);
			M2_tot += pow(M,2.);
			M_tot += M/double(num_rows*num_rows);			
		}

		// Check magnetization stabalization
		// if(i%500 == 0){
		// 	std::cout << num_rows << "\t" << temp << "\t" << i << "\t" << M_tot/current_step << std::endl;
		// }
	}
}

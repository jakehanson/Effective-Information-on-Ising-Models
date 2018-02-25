#include "header.h"

/* Function to calculate the Effective Information */
double get_EI(std::vector<double> ID,std::vector<double> ED,std::vector<std::vector<double>> TPM){
	
	double sum = 0.;
	for(int i = 0;i<TPM.size();i++){
		double D_KL = 0.;
		for(int j = 0;j<TPM.size();j++){
			// Only count non-zero TP
			if(TPM[i][j] > 0.){
				D_KL += TPM[i][j]*std::log2(TPM[i][j]/ED[j]);
			}
		}
		sum = sum + D_KL*ID[i];
	}

	return sum;
}
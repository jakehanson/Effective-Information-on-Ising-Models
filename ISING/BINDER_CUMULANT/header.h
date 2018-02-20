#pragma once
#include <vector>
#include <iostream>
#include <random>
#include <stdio.h>
#include <map>


/* Structure to hold the current properties of the Ising Model */
struct Ising_Model
{
	Ising_Model(int n_rows, int n_cols, double T, long N_steps,long N_burns);  // signature for constructor, construct instance of class

	std::vector<std::vector<int>> spin_matrix;	// 2d array for the current spin matrix
	long current_step; // integer denoting the number of steps taken
	int num_rows; // number of rows
	int num_cols; //number of columns
	long num_steps; // number of simulation steps
	long num_burns; // length of burn in period
	double temp; // temperature of the model
	double M_tot;
	double M2_tot;
	double M4_tot;
	
	/* Define Methods */
	void evolve(void);
};
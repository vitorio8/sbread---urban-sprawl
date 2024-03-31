#ifndef INCLUDE_GUARD_GLOBAL_VARIABLES_H
#define INCLUDE_GUARD_GLOBAL_VARIABLES_H
/*
This is to declare global variables. (Definition is done in global_variables.cpp
*/

#include<vector>
#include<Rcpp.h>
#include"transmission_adj_list.h"
using namespace std;
using namespace Rcpp;

/*
N: number of cells
K: number of possible languages
T: number of timesteps (time 0 -> initial known map.  time T-1 -> known map. time 1...T-2 -> unknown)
*/

extern int N, K, T;

extern double mutation_rate;

// Weighted Adjacency List
extern vector<Transmission_Adjacency_List> cell_adj_lists;

//2-D matrix whose [t][i] element shows the state in cell i at time t. -1 if the state is unknown.
extern vector<vector<int>> known_state_timeline;

//2-D matrix whose [t][i] element shows whether the known disturbance occurred in cell i at time t.
//turbulance event is excluded from the computation of posterior probability.
extern vector<vector<bool>> known_disturbance_timeline;

// Sampled histories. Result of MCMC
extern vector<NumericMatrix> sampled_state_distribution_timelines;


/*
Global variables related to the configuration of MCMC.
Num_iteration: total number of iterations for each chain.
*/

extern long long Num_chains, Length_burn_in, Sample_size, Length_sampling_interval, Num_iterations;


#endif
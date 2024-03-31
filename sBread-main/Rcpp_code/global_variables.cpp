#include<vector>
#include<fstream>
#include"global_variables.h"

using namespace Rcpp;
using namespace std;

//definition of global variables.
//N: number of cells with a language. K: number of possible languages. T: number of timesteps including first and last ones.
int N, K, T;

double mutation_rate;

vector<Transmission_Adjacency_List> cell_adj_lists;

vector<vector<int>> known_state_timeline;
vector<vector<bool>> known_disturbance_timeline;

vector<NumericMatrix> sampled_state_distribution_timelines;

long long Num_chains, Length_burn_in, Sample_size, Length_sampling_interval, Num_iterations;
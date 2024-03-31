//This is the main routine.

#include<iostream>
#include<vector>
#include<string>
#include<thread>
#include<random>
#include<sys/stat.h>
#include<Rcpp.h>
#include"global_variables.h"
#include"data_interface.h"
#include"transmission_adj_list.h"
#include"MH_algorithm.h"
#include"error_check.h"

using namespace std;
using namespace Rcpp;


static void run_mcmc(MCMC mcmc, int chain, int seed) {
  mcmc.MCMC_main(chain, seed);
}

int main(){
    return 0;
}

// [[Rcpp::export]]
vector<NumericMatrix> sBread(int seed_, int K_, double mutation_rate_, int Num_chains_, int Length_burn_in_, int Sample_size_, int Length_sampling_interval_, 
    List adj_list_, List weight_list_, NumericMatrix known_state_timeline_, bool in_between_disturbance_ = false){
    
    //error check.
    error_check(K_, mutation_rate_, Num_chains_, Length_burn_in_, Sample_size_, Length_sampling_interval_, adj_list_, weight_list_, known_state_timeline_);

    int seed = seed_;

    T = known_state_timeline_.nrow();
    N = known_state_timeline_.ncol();
    K = K_;
    mutation_rate = mutation_rate_;
    Num_chains = Num_chains_;
    Length_burn_in = Length_burn_in_;
    Sample_size = Sample_size_;
    Length_sampling_interval = Length_sampling_interval_;

    Num_iterations = Length_burn_in + Sample_size * Length_sampling_interval;

    //initialize the global variable cell_adj_lists
    cell_adj_lists = construct_adj_lists(Rlist_to_2D_vector<int>(adj_list_), Rlist_to_2D_vector<double>(weight_list_));

    known_state_timeline = Rmatrix_to_2D_vector<int>(known_state_timeline_);

    known_disturbance_timeline = vector<vector<bool>>(T, vector<bool>(N, false));
    if(in_between_disturbance_){
        for(int t = 1; t < T - 1; t++){
            for(int i = 0; i < N; i++){
                //all the known history in-between is considered disturbance.
                known_disturbance_timeline.at(t).at(i) = (known_state_timeline.at(t).at(i) != -1);
            }
        }
    }

    //initialize the vector containing the sampled histories.
    sampled_state_distribution_timelines.resize(Num_chains * Sample_size);

    vector<MCMC> mcmc_engines(Num_chains);
    vector<thread> MCMC_threads;

    mt19937 random(seed);
    uniform_int_distribution<> rand_int(0, pow(10, 9));

    for(int i = 0; i < Num_chains; i++){
        MCMC_threads.emplace_back(
                thread(run_mcmc, mcmc_engines.at(i), i, rand_int(random))
        );
    }
    for(int i = 0; i < Num_chains; i++){
        MCMC_threads.at(i).join();
    }

    return(sampled_state_distribution_timelines);
}
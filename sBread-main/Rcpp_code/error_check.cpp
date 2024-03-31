/*
This is the function to check if the input from the R console is correct.
This function should be called right after the main function is called.
*/

#include<Rcpp.h>
#include"error_check.h"

using namespace Rcpp;

void error_check(int K_, double mutation_rate_, int Num_chains_, int Length_burn_in_, int Sample_size_, 
    int Length_sampling_interval_, List adj_list_, List weight_list_, NumericMatrix known_state_timeline_){

    //errors concerning the range of arguments.
    if(K_ < 1){
        stop("K must be a positive integer.");
    }
    if(Num_chains_ < 1){
        stop("Num_chain must be a positive integer.");
    }
    if(Length_burn_in_ < 1){
        stop("Length_burn_in must be a positive integer.");
    }
    if(Sample_size_ < 1){
        stop("Sample_size must be a positive integer.");
    }
    if(Length_sampling_interval_ < 1){
        stop("Length_sampling_interval must be a positive integer.");
    }
    if(mutation_rate_ < 0.0 || mutation_rate_ > 1.0){
        stop("Mutation rate must be positive and must not exceed 1.0.");
    }

    //errors concerning data sizes.
    int adj_list_length = adj_list_.length();
    int weight_list_length = weight_list_.length();
    int known_state_timeline_ncol = known_state_timeline_.ncol();
    if(adj_list_length != weight_list_length || adj_list_length != known_state_timeline_ncol){
        stop("The number of nodes must be consistent in the three arguments, adj_list, weight_list and known_state_timeline.");
    }

    //errors in wrong input in the known history
    for(int i = 0; i < known_state_timeline_.nrow(); i++){
        for(int j = 0; j < known_state_timeline_.ncol(); j++){
            int state_value = known_state_timeline_(i, j);
            if(state_value < -1 || state_value > K_ - 1){
                stop("State values in 'known_state_timeline' must be either -1, 0, 1,...,K-1.");
            }
        }
    }
}
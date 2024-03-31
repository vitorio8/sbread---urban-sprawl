#ifndef INCLUDE_GUARD_ERROR_CHECK_H
#define INCLUDE_GUARD_ERROR_CHECK_H

#include<Rcpp.h>

using namespace Rcpp;

void error_check(int K_, double mutation_rate_, int Num_chains_, int Length_burn_in_, int Sample_size_, 
    int Length_sampling_interval_, List adj_list_, List weight_list_, NumericMatrix known_state_timeline_);

#endif
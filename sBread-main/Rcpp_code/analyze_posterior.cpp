#include<Rcpp.h>
#include<vector>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
vector<NumericMatrix> compute_posterior(int K, int T, int N, List sampled_histories){
    int num_samples = sampled_histories.size();

    vector<NumericMatrix> result(K);
    for(int k = 0; k < K; k++){
        NumericMatrix mat_tmp(T, N);
        result.at(k) = mat_tmp;
    }

    for(int i = 0; i < num_samples; i++){
        NumericMatrix current_history = sampled_histories[i];
        for(int t = 0; t < T; t++){
            for(int j = 0; j < N; j++){
                int state_id = current_history(t, j);
                result.at(state_id)(t, j) += 1.0 / num_samples;
            }
        }
    }

    return(result);
}
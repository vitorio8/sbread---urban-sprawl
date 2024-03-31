#ifndef INCLUDE_GUARD_MH_ALGORITHM_H
#define INCLUDE_GUARD_MH_ALGORITHM_H

#include<vector>
#include <set>
#include <unordered_map>
using namespace std;

/*
We create the class "MCMC" to gather functions concerning the MH-algorithm.
*/

struct Change {
/*
A struct for storing changes to the state_timeline.
*/
    int time;
    int cell;
    int old_value;
    int new_value;
};

class MCMC{
/*
state_timeline[t][i] ... 2D vector. state ID at time t at cell i.
forward_probabilities[t][i][k] ... 3D vector. probability that cell i learns the state k at time t, given map of time t-1.
*/
    int chain_id;
    int num_accept;
    int n_tries = 100000;

    vector<vector<int>> state_timeline;
    vector<vector<vector<double>>> forward_probabilities;
    vector<vector<vector<double>>> backward_probabilities;
    double log_posterior;

    void initialize_state_timeline(int seed);
    void initialize_forward_probabilities();
    void initialize_backward_probabilities();
    void initialize_log_posterior();

    void mutate_state_timeline(int change_time, int change_cell, int new_state_id);

    double propose_pointwise(mt19937& random, vector<Change>& changes);
    Change randomly_expand_changes(mt19937& random, const set<int>& changed_cells, int time);
    double propose_batch(mt19937& random, vector<Change>& changes, unordered_map<int, vector<double>>& changed_forward_probs_t);
    double calculate_posterior_increment(vector<Change>& changes, unordered_map<int, vector<double>>& changed_forward_probs_t);
    void MH_algorithm(int seed);

    public:
        //main function to conduct MCMC.
        void MCMC_main(int chain, int seed);
};


//auxiliary functions
bool is_known(int t, int i);
bool is_disturbance(int t, int i);

#endif
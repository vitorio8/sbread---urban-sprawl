/*
This is the main function for MCMC.
*/

#include<Rcpp.h>
#include<math.h>
#include<vector>
#include<iostream>
#include<fstream>
#include<random>
#include<set>
#include<unordered_map>

#include"global_variables.h"
#include"MH_algorithm.h"

using namespace Rcpp;
using namespace std;

//random value genertors.
uniform_real_distribution<> rand_real(0, 1);
uniform_int_distribution<> rand_int(0, pow(10, 9));


class SBreadException : public std::exception
{
public:
    SBreadException(const char* err) : std::exception() {}
};


void print_change(Change& change) {
    Rcout << "Change (" << change.time << ", " << change.cell << ") from " << change.old_value << " to " << change.new_value << "\n";
}


void MCMC::initialize_state_timeline(int seed){
/*
This function initializes the 2D-vector whose [t][i] element represents the state at time t in cell i.
*/
    state_timeline = known_state_timeline;

    mt19937 random(seed);
    //fill unknown states to initialize the distribution timeline.
    for(int t = 0; t < T; t++){
        for(int i = 0; i < N; i++){
            if(!is_known(t, i)){
                int state_id = rand_int(random) % K;
                state_timeline.at(t).at(i) = state_id;
            }
        }
    }
}

void MCMC::initialize_forward_probabilities(){
/*
This function computes the 3D-vector whose [t][i][k]-th component represents.
the probability that cell i learns the state k at time t, given map of time t-1.
(given that the disturbance does not happen.)
This function is used at the first iteration of MCMC.

Input: state_timeline[t][i] ... state ID at time t at cell i. 
*/
    vector<vector<vector<double>>> initial_forward_probabilities(T, vector<vector<double>>(N, vector<double>(K, 0.0)));
    forward_probabilities = initial_forward_probabilities;

    for(int i = 0; i < N; i++){
        for(int k = 0; k < K; k++){
            forward_probabilities.at(0).at(i).at(k) = -1.0;
        }
    }

    for(int t = 1; t < T; t++){
        for(int i = 0; i < N; i++){
            //learning the state by transmission.
            for(int j = 0; j < cell_adj_lists.at(i).num_neighbors; j++){
                int model_cell_id = cell_adj_lists.at(i).adjacency_list.at(j);
                int model_state_id = state_timeline.at(t-1).at(model_cell_id);
                forward_probabilities.at(t).at(i).at(model_state_id) += (1 - mutation_rate) * cell_adj_lists.at(i).transmission_rate.at(j);
            }

            //learning the state by mutation.
            for(int k = 0; k < K; k++){
                forward_probabilities.at(t).at(i).at(k) += mutation_rate / K;
            }
        }
    }
}


void MCMC::initialize_backward_probabilities(){
/*
This function computes the 3D-vector whose [t][i][k]-th component represents
the backwards convolutions which are used for the MCMC proposals at each cell i
to have state k at time t, given map of time t+1.
(given that the disturbance does not happen.)
This function is used at the first iteration of MCMC.

Input: state_timeline[t][i] ... state ID at time t at cell i.
*/
    vector<vector<vector<double>>> initial_backward_probabilities(T, vector<vector<double>>(N, vector<double>(K, 0.0)));
    backward_probabilities = initial_backward_probabilities;

    for(int i = 0; i < N; i++){
        for(int k = 0; k < K; k++){
            backward_probabilities.at(T-1).at(i).at(k) = -1.0;
        }
    }

    for(int t = T-2; t >= 0; t--){
        for(int i = 0; i < N; i++){
            //learning the state by transmission.
            for(int j = 0; j < cell_adj_lists.at(i).num_neighbors; j++){
                int model_cell_id = cell_adj_lists.at(i).adjacency_list.at(j);
                int model_state_id = state_timeline.at(t+1).at(model_cell_id);
                backward_probabilities.at(t).at(i).at(model_state_id) += (1 - mutation_rate) * cell_adj_lists.at(i).transmission_rate.at(j);
            }

            //learning the state by mutation.
            for(int k = 0; k < K; k++){
                backward_probabilities.at(t).at(i).at(k) += mutation_rate / K;
            }
        }
    }
}


void MCMC::initialize_log_posterior(){
/*
This function computes the posterior probability expressed as the logarithm with base 10.
Used at the first iteration of MCMC.
*/
    log_posterior = 0.0;
    for(int t = 1; t < T; t++){
        for(int i = 0; i < N; i++){
            if(is_disturbance(t, i)){  //if disturbance, exclude from the computation of posterior.
                continue;
            }

            int state_id = state_timeline.at(t).at(i);
            log_posterior += log10(forward_probabilities.at(t).at(i).at(state_id));
        }
    }
}


void MCMC::mutate_state_timeline(int change_time, int change_cell, int new_state_id){
/*
This function generates a single mutation in the state timeline.
Used when accepting a change in the MH-algorithm.
Change of the learn probability timeline follows.
*/

    int old_state_id = state_timeline.at(change_time).at(change_cell);

    state_timeline.at(change_time).at(change_cell) = new_state_id;

    for(int neighbor_cell_id : cell_adj_lists.at(change_cell).adjacency_list){
        //probability that the neighboring cell copies from the focal cell which changes the state.
        double transmission_rate = cell_adj_lists.at(neighbor_cell_id).neighbor_id_to_transmission_rate[change_cell];
        
        forward_probabilities.at(change_time + 1).at(neighbor_cell_id).at(new_state_id) += (1 - mutation_rate) * transmission_rate;
        forward_probabilities.at(change_time + 1).at(neighbor_cell_id).at(old_state_id) -= (1 - mutation_rate) * transmission_rate;

        backward_probabilities.at(change_time - 1).at(neighbor_cell_id).at(new_state_id) += (1 - mutation_rate) * transmission_rate;
        backward_probabilities.at(change_time - 1).at(neighbor_cell_id).at(old_state_id) -= (1 - mutation_rate) * transmission_rate;
    }
}

double MCMC::propose_pointwise(mt19937& random, vector<Change>& changes) {

    // int change_time, change_cell;
    Change change;
    while(true){
        change.time = rand_int(random) % (T - 2) + 1;  // uniform integer distribution within the range of [1, T-2]
        change.cell = rand_int(random) % N;
        if(!is_known(change.time, change.cell)){  //repeat until the unknown timestep and cell are chosen.
            break;
        }
    }
    change.old_value = state_timeline.at(change.time).at(change.cell);

    while(true){
        change.new_value = rand_int(random) % K;
        if(change.new_value != change.old_value){
            break;
        }
    }

    //Increment of posterior probability measured in logarithmic function of base 10.
    double log_posterior_increment = 0.0;

    //Compute the change of probability to learn a state given the previous timestep
    log_posterior_increment -= log10(forward_probabilities.at(change.time).at(change.cell).at(change.old_value));
    log_posterior_increment += log10(forward_probabilities.at(change.time).at(change.cell).at(change.new_value));

    for(int i = 0; i < cell_adj_lists.at(change.cell).num_neighbors; i++){
        int neighbor_id = cell_adj_lists.at(change.cell).adjacency_list.at(i);

        //probability that the neighboring cell copies from the focal cell which changes the state.
        double transmission_rate = cell_adj_lists.at(neighbor_id).neighbor_id_to_transmission_rate[change.cell];

        int neighbor_state_id = state_timeline.at(change.time + 1).at(neighbor_id);

        if(is_disturbance(change.time, neighbor_state_id)){  //if disturbance, exclude from the computation of posterior.
            continue;
        }

        if(neighbor_state_id == change.old_value){
            log_posterior_increment -= log10(forward_probabilities.at(change.time + 1).at(neighbor_id).at(neighbor_state_id));
            log_posterior_increment += log10(forward_probabilities.at(change.time + 1).at(neighbor_id).at(neighbor_state_id) - (1 - mutation_rate) * transmission_rate);

        }else if(neighbor_state_id == change.new_value){
            log_posterior_increment -= log10(forward_probabilities.at(change.time + 1).at(neighbor_id).at(neighbor_state_id));
            log_posterior_increment += log10(forward_probabilities.at(change.time + 1).at(neighbor_id).at(neighbor_state_id) + (1 - mutation_rate) * transmission_rate);
        }
    }

    changes.push_back(change);

    return log_posterior_increment;
}


int random_neighbor(int node, mt19937& random) {
    Transmission_Adjacency_List neighbors = cell_adj_lists.at(node);
    int i = rand_int(random) % neighbors.num_neighbors;
    return neighbors.adjacency_list.at(i);
}

int random_draw_from_set(const set<int>& s, mt19937& random) {
    double r = rand_int(random) % s.size();
    std::set<int>::iterator it = s.begin();
    for (; r != 0; r--) it++;
    return *it;
}

Change MCMC::randomly_expand_changes(mt19937& random, const set<int>& changed_cells, int time) {
    for (int i_try=0; i_try < n_tries; i_try++) {
        int ref_cell = random_draw_from_set(changed_cells, random);
        int cell = random_neighbor(ref_cell, random);
        if (!is_known(time, cell) && (changed_cells.count(cell) == 0)) {
            return Change {time, cell, -1, -1};
        }
    }
    //for debug
    for(int i = 0; i < 10000; i++){
        Rcout << i << endl;
    }
    throw SBreadException("No neighbours found");
}

int sample_discrete(const std::vector<double>& probabilities, mt19937& random) {
    std::discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());
    return distribution(random);
}

vector<double> multiply(const vector<double>& x, const vector<double>& y) {
    vector<double> xy;
    for (int i=0; i < x.size(); ++i)
        xy.push_back(x[i]*y[i]);
    return xy;
}

double sum(vector<double>& x) {
    double sum_x = 0;
    for (int i=0; i < x.size(); ++i)
        sum_x += x[i];
    return sum_x;
}

vector<double> normalize(vector<double> x) {
    vector<double> xy;
    double sum_x = sum(x);
    for (int i=0; i < x.size(); ++i)
        xy.push_back(x[i] / sum_x);
    return xy ;
}

double MCMC::propose_batch(mt19937& random, vector<Change>& changes, unordered_map<int, vector<double>>& changed_forward_probs_t) {
    // Ratio of the forward and backward probability density of the proposed step change in base 10 log-space
    double log_proposal_prob_ratio = 0.0;

    // int change_time, change_cell;
    Change first_change;
    while(true){
        first_change.time = rand_int(random) % (T - 2) + 1;  // uniform integer distribution within the range of [1, T-2]
        first_change.cell = rand_int(random) % N;
        if(!is_known(first_change.time, first_change.cell)){  //repeat until the unknown timestep and cell are chosen.
            break;
        }
    }
    int time = first_change.time;
    first_change.old_value = state_timeline.at(first_change.time).at(first_change.cell);

    //check the number of neighbours without the known state.
    int num_maximum_changes = 0;
    for(int c : cell_adj_lists.at(first_change.cell).adjacency_list){
        if(!is_known(first_change.time, c)){
            num_maximum_changes += 1;
        }
    }
    int n_changes = 1 + (rand_int(random) % 3);  // Random number from 1 to 3
    n_changes = min(n_changes, num_maximum_changes);  //number of changes is up to the number of neighbours without known state.
    set<int> changed_cells;


    // Remember the change (and the changed cell separately for convenience)
    changes.push_back(first_change);
    changed_cells.insert(first_change.cell);

    // Iteratively choose neighboring cells to add to the changes
    for (int i_change = 0; i_change < n_changes-1; i_change++) {
        Change next_change = randomly_expand_changes(random, changed_cells, time);
        next_change.old_value = state_timeline[time][next_change.cell];
        changes.push_back(next_change);
        changed_cells.insert(next_change.cell);
    }

    // Sample new states for the selected cells at the selected time
    set<int> affected_cells;
    for (auto change = changes.begin(); change != changes.end(); change++) {
        // Create a proposal distribution for the state in (change->cell) at (change->time)
        vector<double> probs_forw = forward_probabilities.at(change->time).at(change->cell);
        vector<double> probs_back = backward_probabilities.at(change->time).at(change->cell);
        vector<double> probs = normalize(multiply(probs_forw, probs_back));
        
        // Sample from the proposal distribution
        change->new_value = sample_discrete(probs, random);
        if (change->new_value >= K)
            throw SBreadException("");

        // Compute the transition probabilities
        log_proposal_prob_ratio -= log10(probs.at(change->new_value));
        log_proposal_prob_ratio += log10(probs.at(change->old_value));

        // Update changed_forward_probs_t to more efficiently caluclate the posterior increment later on
        for (int c : cell_adj_lists.at(change->cell).adjacency_list) {
            // Move on if the state at cell c is not affected by the change
            int c_state = state_timeline.at(time + 1).at(c);
            if (c_state != change->old_value && c_state != change->new_value) {
                continue;
            }

            // Remember c as an affected cell
            affected_cells.insert(c);

            // Register a changed forward_probability for cell c at time t+1
            if (changed_forward_probs_t.count(c) == 0) {
                vector<double> old_forward_prob = forward_probabilities.at(time+1).at(c);
                changed_forward_probs_t.insert({c, old_forward_prob});
            }

            // Calculate the changed forward probability
            double transmission_rate = cell_adj_lists.at(c).neighbor_id_to_transmission_rate[change->cell];
            changed_forward_probs_t[c][change->new_value] += (1 - mutation_rate) * transmission_rate;
            changed_forward_probs_t[c][change->old_value] -= (1 - mutation_rate) * transmission_rate;
        }
    }

    return log_proposal_prob_ratio;
}


double MCMC::calculate_posterior_increment(vector<Change>& changes, unordered_map<int, vector<double>>& changed_forward_probs_t) {
    double log_posterior_increment = 0.0;
    int time = -1;
    //Compute the change of probability to learn a state given the previous timestep
    for (Change change : changes) {
        vector<double> forward_t_c  = forward_probabilities.at(change.time).at(change.cell);
        log_posterior_increment -= log10(forward_t_c.at(change.old_value));
        log_posterior_increment += log10(forward_t_c.at(change.new_value));
        time = change.time;
    }

    for (auto changed_prob : changed_forward_probs_t) {
        int c = changed_prob.first;
        vector<double> new_forward_prob = changed_prob.second;
        int c_state = state_timeline.at(time + 1).at(c);

        if(is_disturbance(time + 1, c)){  //if disturbance, exclude from the computation of posterior.
            continue;
        }
        log_posterior_increment -= log10(forward_probabilities.at(time + 1).at(c).at(c_state));
        log_posterior_increment += log10(new_forward_prob.at(c_state));
    }

    return log_posterior_increment;
}


void MCMC::MH_algorithm(int seed){
/*
This function proposes a new state distribution timeline and compute the local change in the conditional probability.
The proposed change is accepted based on the MH_algorithm.
*/
    vector<Change> changes;
    unordered_map<int, vector<double>> changed_forward_probs_t;
    mt19937 random(seed);

    double log_hgf = propose_batch(random, changes, changed_forward_probs_t);
    double log_posterior_increment = calculate_posterior_increment(changes, changed_forward_probs_t);
    double p_accept = pow(10, log_hgf + log_posterior_increment);

    //Decide whether to accept the proposal.
    if (rand_real(random) < p_accept) {
        for (Change change : changes) {
            mutate_state_timeline(change.time, change.cell, change.new_value);
        }
        log_posterior += log_posterior_increment;
        num_accept += 1;
    }
}


void MCMC::MCMC_main(int chain, int seed){
    chain_id = chain;

    //Random value generators.
    mt19937 random(seed);

    //Initialization.
    initialize_state_timeline(rand_int(random));
    initialize_forward_probabilities();
    initialize_backward_probabilities();
    initialize_log_posterior();

    num_accept = 0;

    int sample_id = 0;  //sample ID for the chain.
    for(long long i = 0; i < Num_iterations; i++){
        if(i % 10000 == 0){
            if(chain_id == 0){
                Rcout << "Iteration: " << i << "   log posterior: " << log_posterior << "   num accepted: " << num_accept << "\n";
            }

            num_accept = 0;
        }

        if(i == Length_burn_in - 1){
            if(chain_id == 0){
                Rcout << "Burn-in is over." << "\n";
            }

            initialize_forward_probabilities();
            initialize_backward_probabilities();
            initialize_log_posterior();
        }

        if(i % Length_sampling_interval == 0 && i >= Length_burn_in){
            if(chain_id == 0){
                //Rcout << "Sampling in progress..." << "\n";
            }

            //sampling the state_timeline.
            NumericMatrix current_sample(T, N);
            for(int t = 0; t < T; t++){
                for(int j = 0; j < N; j++){
                    current_sample(t, j) = state_timeline.at(t).at(j);
                }
            }
            sampled_state_distribution_timelines.at(sample_id + chain_id * Sample_size) = current_sample;
            sample_id += 1;
        }

        MH_algorithm(rand_int(random));
    }
}


//auxiliary functions.

bool is_known(int t, int i){
    //returns true if the state at time t in cell i is known.
    return(known_state_timeline.at(t).at(i) != -1);
}

bool is_disturbance(int t, int i){
    return(known_disturbance_timeline.at(t).at(i));
}
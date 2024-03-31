#include<vector>
#include"transmission_adj_list.h"
#include"global_variables.h"

using namespace std;


vector<Transmission_Adjacency_List> construct_adj_lists(const vector<vector<int>>& adj_list, const vector<vector<double>>& weight_list){
    //Create the vector of Transmission_Adjacency_List class from the adjacency list and weight list.

    vector<Transmission_Adjacency_List> result_adj_lists(N);
    for(int i = 0; i < N; i++){
        result_adj_lists.at(i).adjacency_list = adj_list.at(i);
        result_adj_lists.at(i).transmission_rate = weight_list.at(i);
        result_adj_lists.at(i).num_neighbors = adj_list.at(i).size();

        for(int j = 0; j < result_adj_lists.at(i).num_neighbors; j++){
            result_adj_lists.at(i).neighbor_id_to_transmission_rate[result_adj_lists.at(i).adjacency_list.at(j)] = result_adj_lists.at(i).transmission_rate.at(j);
        }
    }

    return(result_adj_lists);
}
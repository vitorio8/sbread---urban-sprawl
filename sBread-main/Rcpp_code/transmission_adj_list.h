#ifndef INCLUDE_GUARD_TRANSMISSION_ADJ_LIST_H
#define INCLUDE_GUARD_TRANSMISSION_ADJ_LIST_H

/* This is the struct which treats the adjacency list of each population (cell) weighted by transmission rates.
   We need this adjacency list for every cell.
*/

#include<vector>
#include<unordered_map>
using namespace std;

struct Transmission_Adjacency_List{
    int num_neighbors;
    vector<int> adjacency_list;
    vector<double> transmission_rate;
    unordered_map<int, double> neighbor_id_to_transmission_rate;
};

vector<Transmission_Adjacency_List> construct_adj_lists(const vector<vector<int>>& adj_list, const vector<vector<double>>& weight_list);

#endif
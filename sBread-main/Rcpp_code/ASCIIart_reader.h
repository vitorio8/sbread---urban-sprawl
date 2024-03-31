#ifndef INCLUDE_ASCIIART_READER_H
#define INCLUDE_ASCIIART_READER_H
#include<vector>
#include<string>
#include<unordered_map>
using namespace std;

struct Landscape{
    /*
    cell_id[x][y]: 2D-vector. cell ID displayed on the 2D plein.
    is_land[x][y]: True if cell at (x, y) is land, and false if cell at (x, y) is sea.
    */

    vector<vector<int>> cell_id;
    vector<vector<bool>> is_land;

    int x_size, y_size;

    void initialize_landscape(vector<vector<vector<int>>> initial_timeline_txy);
};

struct Transmission_Adjacency_List{
    int num_neighbors;
    vector<int> adjacency_list;
    vector<double> transmission_rate;
    unordered_map<int, double> neighbor_id_to_transmission_rate;
};


//prototype declaration
void initialize_xyN(string input_directory_name);
vector<vector<vector<int>>> initial_timeline_txy_from_input(string input_directory_name);
vector<vector<int>> initial_language_distribution_timeline_from_timeline_txy(const vector<vector<vector<int>>>& initial_timeline_txy);
vector<Transmission_Adjacency_List> construct_adj_lists(const Landscape& landscape);


#endif
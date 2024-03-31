#include<fstream>
#include<iostream>
#include<vector>
#include<string>
#include<Rcpp.h>
#include"ASCIIart_reader.h"
using namespace std;
using namespace Rcpp;

//definition of global variables.
namespace global_variables{
    int N, T;
    int x_size, y_size;

    Landscape landscape_xy;
    vector<Transmission_Adjacency_List> cell_adj_lists;

    int NoLand_Code = -100;
    int State_Unknown_Code = -1;
}
namespace gv = global_variables;


void initialize_xyN(string input_directory_name){
/*
This is to initialize three global variables x_size, y_size, and N.    
*/
    string file_name = input_directory_name + "/T0.txt";
    fstream reading_file(file_name);
    string reading_line_buffer;

    int y_count = 0;
    int N_count = 0;
    while(!reading_file.eof()){
        y_count += 1;
        getline(reading_file, reading_line_buffer);

        for(int i = 0; i < reading_line_buffer.size(); i++){
            if(reading_line_buffer.at(i) != '.'){
                N_count += 1;
            }
        }
    }

    gv::x_size = reading_line_buffer.size();
    gv::y_size = y_count;
    gv::N = N_count;
}


vector<vector<vector<int>>> initial_timeline_txy_from_input(string input_directory_name){
/*
This function returns the initial language_distribution_timeline in the 3D (t-x-y) format.
initial_timeline_txy[t][x][y] : Language at time t at (x, y)
markers:  -1 -> no language   -2 -> masked (known settlement)

It is necessary to initialize the global variables x_size and y_size in advance.
*/

    //result 3D-vector
    vector<vector<vector<int>>> initial_timeline_txy(gv::T, vector<vector<int>>(gv::x_size, vector<int>(gv::y_size, gv::State_Unknown_Code)));

    for(int t = 0; t < gv::T; t++){
        string current_file_name = input_directory_name + "/T" + to_string(t) + ".txt";

        fstream initial_data_file(current_file_name);

        if(!initial_data_file){  //if the input file does not exist, copy the places of no-land. States of all other cells are unknown.
            for(int y = 0; y < gv::y_size; y++){
                for(int x = 0; x < gv::x_size; x++){
                    if(initial_timeline_txy.at(t-1).at(x).at(y) == gv::NoLand_Code){
                        initial_timeline_txy.at(t).at(x).at(y) = gv::NoLand_Code;
                    }
                }
            }
            continue;  //move to the next t.
        }

        string reading_line_buffer;

        for(int y = 0; y < gv::y_size; y++){
            getline(initial_data_file, reading_line_buffer);

            for(int x = 0; x < gv::x_size; x++){
                if(reading_line_buffer.at(x) == '.'){  //not land
                    initial_timeline_txy.at(t).at(x).at(y) = gv::NoLand_Code;
                }else if(reading_line_buffer.at(x) == '?'){  //Unknown
                    initial_timeline_txy.at(t).at(x).at(y) = gv::State_Unknown_Code;
                }else{
                    initial_timeline_txy.at(t).at(x).at(y) = reading_line_buffer.at(x) - '0'; //char -> int conversion
                }
            }
        }
    }

    return(initial_timeline_txy);
}


//Definition of the function in the Landscape struct.
void Landscape::initialize_landscape(vector<vector<vector<int>>> initial_timeline_txy){
/*
Create the Landscape struct from the initial_language_distribution_timiline given in xyt format.
initial_timeline_txy[t][x][y] : language at time t at (x, y).
*/

    x_size = initial_timeline_txy.at(0).size();
    y_size = initial_timeline_txy.at(0).at(0).size();

    cell_id.resize(x_size);
    is_land.resize(x_size);
    for(int x = 0; x < x_size; x++){
        cell_id.at(x).resize(y_size);
        is_land.at(x).resize(y_size);
    }

    int current_cell_id = 0;
    for(int x = 0; x < x_size; x++){
        for(int y = 0; y < y_size; y++){
            if(initial_timeline_txy.at(0).at(x).at(y) != gv::NoLand_Code){
                cell_id.at(x).at(y) = current_cell_id;
                is_land.at(x).at(y) = true;
                current_cell_id += 1;
            }else{
                cell_id.at(x).at(y) = -1;
                is_land.at(x).at(y) = false;                
            }
        }
    }
}


vector<vector<int>> initial_language_distribution_timeline_from_timeline_txy(const vector<vector<vector<int>>>& initial_timeline_txy){
/*
This is to convert the language distribution timeline of 3D format (t-x-y) into the 2D format (t-cell_id).
It is necessary to initialize landscape_xy before calling this function.
*/
    vector<vector<int>> initial_language_distribution_timeline(gv::T, vector<int>(gv::N, gv::State_Unknown_Code));

    for(int x = 0; x < gv::x_size; x++){
        for(int y = 0; y < gv::y_size; y++){
            if(gv::landscape_xy.is_land.at(x).at(y)){
                int cell_id = gv::landscape_xy.cell_id.at(x).at(y);

                for(int t = 0; t < gv::T; t++){
                    initial_language_distribution_timeline.at(t).at(cell_id) = initial_timeline_txy.at(t).at(x).at(y);
                    if(initial_timeline_txy.at(t).at(x).at(y) == gv::State_Unknown_Code){
                        initial_language_distribution_timeline.at(t).at(cell_id) = gv::State_Unknown_Code;
                    }
                }
            }
        }
    }

    return(initial_language_distribution_timeline);
}


vector<Transmission_Adjacency_List> construct_adj_lists(const Landscape& landscape){
/*
This is to create a vector of transmission_adj_lists. (i.e., weighted adjacency list for every cell.)
Every cell learns from 8 surrounding cells and self with the equal probability (1/9).
*/
    vector<Transmission_Adjacency_List> result_adj_lists(gv::N);

    for(int x = 0; x < gv::x_size; x++){
        for(int y = 0; y < gv::y_size; y++){
            if(!landscape.is_land.at(x).at(y)){  //if sea, skip.
                continue;
            }

            int cell_id = landscape.cell_id.at(x).at(y);
            
            vector<int> current_adj_list;

            //learning from self
            current_adj_list.push_back(cell_id);

            if(y != 0 && landscape.is_land.at(x).at(y-1)){
                int model_cell_id = landscape.cell_id.at(x).at(y-1);
                current_adj_list.push_back(model_cell_id);
            }
            if(y != gv::y_size - 1 && landscape.is_land.at(x).at(y+1)){
                int model_cell_id = landscape.cell_id.at(x).at(y+1);
                current_adj_list.push_back(model_cell_id);
            }
            
            if(x != 0){
                if(y != 0 && landscape.is_land.at(x-1).at(y-1)){
                    int model_cell_id = landscape.cell_id.at(x-1).at(y-1);
                    current_adj_list.push_back(model_cell_id);
                }

                if(landscape.is_land.at(x-1).at(y)){
                    int model_cell_id = landscape.cell_id.at(x-1).at(y);
                    current_adj_list.push_back(model_cell_id);
                }

                if(y != gv::y_size - 1 && landscape.is_land.at(x-1).at(y+1)){
                    int model_cell_id = landscape.cell_id.at(x-1).at(y+1);
                    current_adj_list.push_back(model_cell_id);
                }
            }

            if(x != gv::x_size - 1){
                if(y != 0 && landscape.is_land.at(x+1).at(y-1)){
                    int model_cell_id = landscape.cell_id.at(x+1).at(y-1);
                    current_adj_list.push_back(model_cell_id);
                }

                if(landscape.is_land.at(x+1).at(y)){
                    int model_cell_id = landscape.cell_id.at(x+1).at(y);
                    current_adj_list.push_back(model_cell_id);
                }

                if(y != gv::y_size - 1 && landscape.is_land.at(x+1).at(y+1)){
                    int model_cell_id = landscape.cell_id.at(x+1).at(y+1);
                    current_adj_list.push_back(model_cell_id);
                }
            }

            result_adj_lists.at(cell_id).adjacency_list = current_adj_list;
            result_adj_lists.at(cell_id).num_neighbors = current_adj_list.size();

            vector<double> transmission_rate_vec(current_adj_list.size(), 1.0 / current_adj_list.size());
            result_adj_lists.at(cell_id).transmission_rate = transmission_rate_vec;

            //assign the unordered_map [neighbor_id -> transmission_rate]
            for(int i = 0; i < current_adj_list.size(); i++){
                int neighbor_id = current_adj_list.at(i);
                double transmission_rate = transmission_rate_vec.at(i);
                result_adj_lists.at(cell_id).neighbor_id_to_transmission_rate[neighbor_id] = transmission_rate;
            }
        }
    }

    return(result_adj_lists);
}


// [[Rcpp::export]]
List ASCIIart_reader(string input_directory_name, int T_){
/*
This is the function called by R.
- initialize the values of x_size, y_size, and N.
- create the language distribution timeline in the 3D format (t-x-y).
- construct the global variable landscape_xy.
- construct the global variable mask_list.
- construct the global variable initial_language_distribution_timeline.
*/

    gv::T = T_; //initialize the global variable.
    initialize_xyN(input_directory_name);
    vector<vector<vector<int>>> initial_timeline_txy = initial_timeline_txy_from_input(input_directory_name);

    gv::landscape_xy.initialize_landscape(initial_timeline_txy);

    //after initializing landscape_xy, we initialize the language_distribution_timeline in the 2D (t-cell_id) format.
    vector<vector<int>> initial_language_distribution_timeline = initial_language_distribution_timeline_from_timeline_txy(initial_timeline_txy);

    //construct the weighted adjacency list.
    gv::cell_adj_lists = construct_adj_lists(gv::landscape_xy);

    //Here we construct Rcpp objects.
    //*
    NumericMatrix known_data(gv::T, gv::N);
    for(int t = 0; t < gv::T; t++){
        for(int i = 0; i < gv::N; i++){
            known_data(t, i) = initial_language_distribution_timeline.at(t).at(i);
        }
    }

    vector<vector<int>> adj_list(gv::N);
    vector<vector<double>> weight_list(gv::N);
    for(int i = 0; i < gv::N; i++){
        adj_list.at(i) = gv::cell_adj_lists.at(i).adjacency_list;
        weight_list.at(i) = gv::cell_adj_lists.at(i).transmission_rate;
    }

    List result_list = List::create(known_data, adj_list, weight_list);
    return(result_list);
    //*/
}

int main(){
    ASCIIart_reader("test", 2);
}
#ifndef INCLUDE_DATA_INTERFACE_VARIABLES_H
#define INCLUDE_DATA_INTERFACE_VARIABLES_H
using namespace std;
using namespace Rcpp;

/*
Functions to convert Rcpp objects into the std objects.
*/

template<typename T1>
vector<vector<T1>> Rlist_to_2D_vector(const List& r_list);

template<typename T1>
vector<vector<T1>> Rmatrix_to_2D_vector(const NumericMatrix& r_mat);

#endif
#include<Rcpp.h>
#include<vector>
#include"data_interface.h"
using namespace std;
using namespace Rcpp;

template<typename T1>
vector<vector<T1>> Rlist_to_2D_vector(const List& r_list){
    T1 N = r_list.size();
    vector<vector<T1>> result_mat(N);

    for(int i = 0; i < N; i++){
        NumericVector tmp_vec = r_list[i];
        int L = tmp_vec.size();
        for(int j = 0; j < L; j++){
            result_mat.at(i).push_back(tmp_vec[j]);
        }
    }

    return(result_mat);
}
//explicit instantiation follows.
template vector<vector<int>> Rlist_to_2D_vector(const List& r_list);
template vector<vector<double>> Rlist_to_2D_vector(const List& r_list);


template<typename T1>
vector<vector<T1>> Rmatrix_to_2D_vector(const NumericMatrix& r_mat){
    int N = r_mat.nrow();
    int M = r_mat.ncol();
    vector<vector<T1>> result_mat(N, vector<T1>(M));

    for(int i = 0; i < N; i++){
        for(int j = 0; j < M; j++){
            result_mat.at(i).at(j) = r_mat(i, j);
        }
    }

    return(result_mat);
}
//explicit instantiation follows.
template vector<vector<int>> Rmatrix_to_2D_vector(const NumericMatrix& r_mat);
template vector<vector<double>> Rmatrix_to_2D_vector(const NumericMatrix& r_mat);
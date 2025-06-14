#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
MatrixXd edge_dist(MatrixXi edge_index,
                   MatrixXd Y) {
  // Initialization
  int n_choose_2(edge_index.rows());
  int i(0);
  int j(0);
  
  MatrixXd edge_dist(n_choose_2, 1);
  
  // Main
  for (int edge_n = 0; edge_n < n_choose_2; ++edge_n) {
    i = edge_index(edge_n, 0) - 1;
    j = edge_index(edge_n, 1) - 1;
    
    edge_dist(edge_n, 0) = abs(Y(i, 0) - Y(j, 0));
  }
  
  return edge_dist;
}

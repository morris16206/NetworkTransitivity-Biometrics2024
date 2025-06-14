#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
MatrixXd alpha_estimation(MatrixXd r) {
  // Initialization
  int n_choose_k(r.rows());
  int p(r.cols());
  MatrixXd alpha_mat(MatrixXd(p, p).setZero());
  
  // Main
  for (int i = 0; i < n_choose_k; ++i) {
    for (int j = 0; j < p; ++j) {
      for (int k = 0; k < p; ++k) {
        if (j != k) {
          alpha_mat(j, k) += r(i, j) * r(i, k);
        }
      }
    }
  }
  
  return alpha_mat / n_choose_k;
}

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
SEXP Gauss_newton_triangle(VectorXd mr1,
                           VectorXd mr2,
                           MatrixXd mw,
                           MatrixXd mm,
                           MatrixXd beta1,
                           MatrixXd beta2,
                           VectorXd pi,
                           int n_time,
                           MatrixXd R) {
  // Initialization
  int p(mm.cols());
  int n_choose_k(mr1.size());
  int restart(0);
  
  VectorXd extended_pi(VectorXd(n_time).setOnes());
  VectorXd square_pi(n_time-1);
  VectorXd cube_pi(n_time);
  VectorXd w(VectorXd(n_time).setOnes());
  VectorXd h(3*n_time-1);
  VectorXd X(p);
  VectorXd f1_tmp(n_time);
  VectorXd f2_tmp(n_time);
  VectorXd f(3*n_time-1);
  MatrixXd D(MatrixXd(3*n_time-1, 3*n_time-1).setIdentity());
  MatrixXd A_sqrt(MatrixXd(3*n_time-1, 3*n_time-1).setZero());
  VectorXd S(3*n_time-1);
  MatrixXd V(3*n_time-1, 3*n_time-1);
  MatrixXd V_inv(3*n_time-1, 3*n_time-1);
  VectorXd U(VectorXd(3*n_time-1).setZero());
  MatrixXd dU(MatrixXd(3*n_time-1, 3*n_time-1).setZero());
  
  // Main
  if (n_time != 1) {
    if ((pi.array() < 1.5e-08 || pi.array() > 1 - 1.5e-08).any()) {
      goto no_update;
    }
    extended_pi.head(n_time-1) = pi;
    square_pi = pi.array().square();
    cube_pi = extended_pi.array().cube();
    f.tail(n_time-1) = pi;
    A_sqrt.diagonal().tail(n_time-1) = (pi.array() * (1 - pi.array())).sqrt();
  }
  
  for (int comb_n = 0; comb_n < n_choose_k; ++comb_n) {
    if (n_time != 1) {
      w.head(n_time-1) = mw.row(comb_n).array().floor();
      h.head(n_time) = w * mr1[comb_n];
      h.segment(n_time, n_time) = w * mr2[comb_n];
      h.tail(n_time-1) = mw.row(comb_n);
    } else {
      h(0) = mr1[comb_n];
      h(1) = mr2[comb_n];
    }
    X = mm.row(comb_n);
    ArrayWrapper Xbeta1_exp((-beta1 * X).array().exp());
    ArrayWrapper Xbeta2_exp((-beta2 * X).array().exp());
    f1_tmp = 1 / (1 + Xbeta1_exp);
    f2_tmp = 1 / (1 + Xbeta2_exp);
    if ((f1_tmp.array() < 1.5e-08 || f1_tmp.array() > 1 - 1.5e-08 || f2_tmp.array() < 1.5e-08 || f2_tmp.array() > 1 - 1.5e-08).any()) {
      goto no_update;
    }
    if (n_time != 1) {
      f.head(n_time) = cube_pi.array() * f1_tmp.array() * f2_tmp.array();
      f.segment(n_time, n_time) = cube_pi.array() * f2_tmp.array();
      D.diagonal().head(n_time) = cube_pi.array() * ((Xbeta1_exp / (1 + Xbeta1_exp).square()).matrix() * X).array() * f2_tmp.array();
      D.diagonal().segment(n_time, n_time) = cube_pi.array() * ((Xbeta2_exp / (1 + Xbeta2_exp).square()).matrix() * X).array();
      D.diagonal(n_time).segment(0, n_time) = cube_pi.array() * ((Xbeta2_exp / (1 + Xbeta2_exp).square()).matrix() * X).array() * f1_tmp.array();
      D.diagonal(n_time).segment(n_time, n_time-1) = 3 * square_pi.array() * f2_tmp.segment(0, n_time-1).array();
      D.diagonal(2*n_time) = 3 * square_pi.array() * f1_tmp.segment(0, n_time-1).array() * f2_tmp.segment(0, n_time-1).array();
      A_sqrt.diagonal().head(n_time) = cube_pi.array() * (f1_tmp.array() * f2_tmp.array() * (1 - f1_tmp.array() * f2_tmp.array())).sqrt();
      A_sqrt.diagonal().segment(n_time, n_time) = cube_pi.array() * (f2_tmp.array() * (1 - f2_tmp.array())).sqrt();
    } else {
      f(0) = (f1_tmp * f2_tmp)[0];
      f(1) = f2_tmp[0];
      D(0, 0) = ((Xbeta1_exp / (1 + Xbeta1_exp).square()).matrix() * X * f2_tmp)[0];
      D(0, 1) = ((Xbeta2_exp / (1 + Xbeta2_exp).square()).matrix() * X * f1_tmp)[0];
      D(1, 1) = ((Xbeta2_exp / (1 + Xbeta2_exp).square()).matrix() * X)[0];
      A_sqrt.diagonal() = (f.array() * (1 - f.array())).sqrt();
    }
    
    S = h - f;
    V = A_sqrt * R * A_sqrt;
    V_inv = V.inverse();
    
    U += D.transpose() * V_inv * S;
    dU += D.transpose() * V_inv * D;
  }
  U /= n_choose_k;
  dU /= n_choose_k;
  
  // Update
  if (restart == 0) {
    MatrixXd update(dU.inverse() * U);
    MatrixXd beta1_update(beta1 + update.topRows(n_time));
    MatrixXd beta2_update(beta2 + update.middleRows(n_time, n_time));
    if (n_time != 1) {
      VectorXd pi_update(pi + update.bottomRows(n_time-1));
      return List::create(Named("beta1") = wrap(beta1_update),
                          Named("beta2") = wrap(beta2_update),
                          Named("pi") = wrap(pi_update));
    } else {
      return List::create(Named("beta1") = wrap(beta1_update),
                          Named("beta2") = wrap(beta2_update),
                          Named("pi") = wrap(pi));
    }
  } else {
    no_update:
    return List::create(Named("beta1") = wrap(beta1),
                        Named("beta2") = wrap(beta2),
                        Named("pi") = wrap(pi));
  }
}

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "do_rgig1.h"
#include <math.h>
using namespace Rcpp;


void res_protector(double& x){
  if (std::abs(x) < DBL_MIN * std::pow(10, 10)){
    double sign = std::copysign(1, x);
    x = DBL_MIN * std::pow(10, 10) * sign;
  }
}


void sample_beta_tilde(arma::mat& beta_nc_samp, arma::vec& y, arma::mat& x, arma::colvec& theta_sr,
                       arma::vec& sig2, arma::colvec& beta_mean, int N, int d, arma::vec delta, Function Rchol){ //
  
  arma::vec y_star = y - x * beta_mean;
  arma::colvec theta = arma::pow(theta_sr, 2);
  arma::mat theta_sr_diag = arma::diagmat(theta_sr);
  arma::mat Ft = x * theta_sr_diag;
  //arma::mat delta_m = arma::diagmat(arma::pow(delta, -0.5));
  arma::mat I_d = arma::eye(d, d);
  
  arma::mat ft;
  arma::mat Qt;
  arma::mat At;
  arma::mat at(d, N+1, arma::fill::zeros);
  arma::mat mt(d, N+1, arma::fill::zeros);
  arma::cube Rt(d, d, N+1, arma::fill::zeros);
  arma::cube Ct(d, d, N+1, arma::fill::zeros);
  arma::mat Bt;
  arma::mat Rtp1_inv;
  Ct.slice(0) = I_d;
  arma::mat mT(d, N+1, arma::fill::zeros);
  arma::cube CT(d, d, N+1, arma::fill::zeros);
  arma::mat L_upper;
  for(int t = 1; t < N+1; t++){
    at.col(t) = mt.col(t-1);
    //Rt.slice(t) = delta_m*Ct.slice(t-1)*delta_m;
    Rt.slice(t) = Ct.slice(t-1) + I_d;
    ft = Ft.row(t-1)*at.col(t);
    Qt = arma::as_scalar(Ft.row(t-1)*Rt.slice(t)*arma::trans(Ft.row(t-1))) + sig2(t-1);
    At = Rt.slice(t)*arma::trans(Ft.row(t-1))/arma::as_scalar(Qt);
    mt.col(t) = at.col(t) + At*(y_star(t-1) - ft);
    Ct.slice(t) = Rt.slice(t) - arma::as_scalar(Qt) * (At*arma::trans(At));
  }
  mT.col(N) = mt.col(N);
  CT.slice(N) = Ct.slice(N);
  //CT.slice(N) = 0.5*CT.slice(N) + 0.5*arma::trans(CT.slice(N));
  arma::mat eps = arma::randn(1, d);
  bool chol_success = chol(L_upper, CT.slice(N));
  // Fall back on Rs chol if armadillo fails
  if(chol_success == false){
    Rcpp::NumericMatrix tmp = Rchol(CT.slice(N), true, false, -1);
    arma::uvec piv = arma::sort_index(as<arma::vec>(tmp.attr("pivot")));
    arma::mat L_upper_tmp = arma::mat(tmp.begin(), d, d, false);
    L_upper = L_upper_tmp.cols(piv);
  }
  //beta_nc_samp.col(N) = mT.col(N) + arma::trans(eps * arma::chol(CT.slice(N)));
  beta_nc_samp.col(N) = mT.col(N) + arma::trans(eps * L_upper);
  for(int t = N-1; t >= 0; t--){
    //Rtp1_inv = arma::inv(Rt.slice(t+1));
    //if(!(Rt.slice(t+1)).is_sympd()){
    //  Rt.slice(t+1) = .5*Rt.slice(t+1) + 0.5*arma::trans(Rt.slice(t+1));
    //}
    Rt.slice(t+1) = .5*Rt.slice(t+1) + 0.5*arma::trans(Rt.slice(t+1));
    Rtp1_inv = arma::inv(Rt.slice(t+1));
    //Rtp1_inv = I_d;
    Bt = Ct.slice(t) * Rtp1_inv;
    //mT.col(t) = mt.col(t) + Bt*(mT.col(t+1) - at.col(t+1));
    mT.col(t) = mt.col(t) + Bt*(beta_nc_samp.col(t+1) - at.col(t+1));
    CT.slice(t) = Ct.slice(t) - Bt*(Rt.slice(t+1))*arma::trans(Bt);
    //CT.slice(t) = 0.5*CT.slice(t) + 0.5*arma::trans(CT.slice(t));
    eps = arma::randn(1, d);
    
    chol_success = chol(L_upper, CT.slice(t));
    if(chol_success == false){
      Rcpp::NumericMatrix tmp = Rchol(CT.slice(N), true, false, -1);
      arma::uvec piv = arma::sort_index(as<arma::vec>(tmp.attr("pivot")));
      arma::mat L_upper_tmp = arma::mat(tmp.begin(), d, d, false);
      L_upper = L_upper_tmp.cols(piv);
    }
    beta_nc_samp.col(t) = mT.col(t) + arma::trans(eps * L_upper);
    //beta_nc_samp.col(t) = mT.col(t);
  }
}


void sample_alpha(arma::vec& alpha_samp, arma::vec& y, arma::mat& x, arma::mat& W, arma::colvec& tau2, arma::colvec& xi2, arma::vec& sig2, arma::vec& a0, int d, Function Rchol) {
  arma::mat W_til = W.t() * arma::diagmat(1/sig2);
  arma::mat A0_sr = arma::diagmat(arma::sqrt(arma::join_cols(tau2, xi2)));
  arma::vec prior_a = a0/arma::join_cols(tau2, xi2);
  arma::mat a = W_til * y + prior_a;
  arma::mat Omega_star = A0_sr * W_til * W * A0_sr + arma::eye(2*d, 2*d);
  
  arma::mat A_t;
  arma::mat A_t_til;
  bool solved = arma::solve(A_t_til, Omega_star, A0_sr, arma::solve_opts::no_approx);
  if (solved == 1){
    A_t = A0_sr * A_t_til;
  } else {
    arma::mat A0 = arma::diagmat(arma::join_cols(tau2, xi2));
    A_t = arma::inv(W_til * W + arma::inv(diagmat(A0)));
  }
  
  arma::vec v = rnorm(2*d);
  
  arma::mat cholA;
  bool chol_success = chol(cholA, A_t);
  
  // Fall back on Rs chol if armadillo fails (it suppports pivoting)
  if (chol_success == false){
    Rcpp::NumericMatrix tmp = Rchol(A_t, true, false, -1);
    arma::mat cholA_tmp = arma::mat(tmp.begin(), 2*d, 2*d, false);
    arma::uvec piv = arma::sort_index(as<arma::vec>(tmp.attr("pivot")));
    cholA = cholA_tmp.cols(piv);
  }
  
  alpha_samp = A_t * a + cholA.t() * v;
  std::for_each(alpha_samp.begin(), alpha_samp.end(), res_protector);
}


void resample_alpha_diff(arma::mat& alpha_samp, arma::mat& betaenter, arma::vec& theta_sr, arma::vec& beta_mean, arma::mat& beta_diff,  arma::vec& xi2, arma::vec& tau2, int d, int N){
  arma::vec sign_sqrt = arma::sign(theta_sr);
  arma::colvec theta_sr_new(d, arma::fill::none);
  int p1_theta = -N/2;
  
  arma::vec theta(d, arma::fill::none);
  
  
  for (int j = 0; j < d; j++){
    double p2_theta = 1/xi2(j);
    double p3_theta = arma::as_scalar(arma::accu(arma::pow(beta_diff.row(j), 2))) +
      std::pow((betaenter(j, 0) - beta_mean(j)), 2);
    
    double res = do_rgig1(p1_theta, p3_theta, p2_theta);
    theta(j) = res;
    theta_sr_new(j) = std::sqrt(arma::as_scalar(theta(j))) * sign_sqrt(j);
  }
  
  arma::colvec beta_mean_new(d, arma::fill::none);
  
  for (int j = 0; j < d; j++){
    double sigma2_beta_mean = 1/(1/tau2(j) + 1/(theta(j)));
    double mu_beta_mean = betaenter(j, 0) * tau2(j)/(tau2(j) + theta(j));
    beta_mean_new(j) = R::rnorm(mu_beta_mean, std::sqrt(sigma2_beta_mean));
  }
  
  alpha_samp = arma::join_cols(beta_mean_new, theta_sr_new);
  std::for_each(alpha_samp.begin(), alpha_samp.end(), res_protector);
}





double sample_kappa2(arma::vec& xi2, double a_xi, double d1, double d2, int d){
  double d1_full = d1 + a_xi * d;
  double d2_full = d2 + arma::as_scalar(arma::mean(xi2)) * a_xi * d * 0.5;
  double kappa2 = R::rgamma(d1_full, 1/d2_full);
  return(kappa2);
}


double sample_lambda2(arma::vec& tau2, double a_tau, double e1, double e2, int d){
  double e1_full = e1 + a_tau * d;
  double e2_full = e2 + arma::as_scalar(arma::mean(tau2)) * a_tau * d * 0.5;
  double lambda2 = R::rgamma(e1_full, 1/e2_full);
  return(lambda2);
}


void sample_sigma2(arma::vec& sig2_samp, arma::vec& y, arma::mat& W, arma::vec& alpha, double c0, double C0, int N){
  double a_full = c0 + N/2;
  double b_full = C0 + 0.5 * arma::as_scalar(arma::sum(arma::pow((y - W*alpha), 2)));
  double sig2 = 1/R::rgamma(a_full, 1/b_full);
  sig2_samp.fill(sig2);
}


double sample_C0(arma::vec& sig2, double g0, double c0, double G0){
  double a_full = g0 + c0;
  double b_full = G0 + 1/sig2(0);
  double C0 = R::rgamma(a_full, 1/b_full);
  return(C0);
}


void sample_eta2(arma::vec& eta2_samp, arma::vec& alpha_samp,
                   arma::vec& delta_samp, double d1, double d2, int d, double par_c){
  double d1_full = d1 + 0.5;
  arma::vec alpha2 = arma::pow(alpha_samp, 2);
  for (int j = 0; j < 2*d; j++){
    double d2_full = d2 + alpha2(j)/(delta_samp(j) + (1 - delta_samp(j)) * par_c);
    double res = 1/R::rgamma(d1_full, 1/d2_full);
    res_protector(res);
    eta2_samp(j) = res;
  }
}


void sample_p(arma::vec& delta_samp, double d1, double d2, int d, double& par_p){
  double delta_sum = arma::sum(delta_samp);
  double d1_full = d1 + delta_sum;
  double d2_full = d2 + 2*d - delta_sum;
  par_p = R::rbeta(d1_full, d2_full);
}


void comp_prob(arma::vec& alpha_samp, arma::vec& eta2_samp, arma::vec& prob_samp, int d, double par_c, double par_p){
  arma::vec eta = arma::pow(eta2_samp, 0.5);
  double par_c_sr = std::sqrt(par_c);
  for (int j = 0; j < 2*d; j++){
    //double tmp1 = par_p * R::dnorm(alpha_samp(j), 0, eta(j), false);
    //double tmp2 = (1 - par_p) * R::dnorm(alpha_samp(j), 0, par_c_sr*eta(j), false);
    //prob_samp(j) = tmp1/(tmp1 + tmp2);
    double l0 = R::dnorm(alpha_samp(j), 0, par_c_sr*eta(j), true) + 1.0e-20;
    double l1 = R::dnorm(alpha_samp(j), 0, eta(j), true) + 1.0e-20;
    prob_samp(j) = 1/(1 + ((1-par_p)/par_p) * exp(l0 - l1));
  }
}


void sample_delta(arma::vec& prob_samp, arma::vec& delta_samp, int d){
  for (int j = 0; j < 2*d; j++){
    delta_samp(j) = R::rbinom(1, prob_samp(j));
  }
}

void comp_SAVS(arma::vec& alpha_samp, arma::mat& W, arma::vec& gamma_samp, int d){
  arma::vec sign_alpha = arma::sign(alpha_samp);
  for(int j = 0; j < 2*d; j++){
    double Z_j = arma::norm(W.col(j), 2);
    double tmp = std::abs(alpha_samp(j))*std::pow(Z_j, 2) - std::pow(std::abs(alpha_samp(j)), -2);
    if (tmp < 0.0){
      tmp = 0.0;
    }
    gamma_samp(j) = sign_alpha(j) * tmp / std::pow(Z_j, 2);
  }
}


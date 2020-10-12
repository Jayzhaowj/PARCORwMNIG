#ifndef SAMPLE_PARAMETERS_H
#define SAMPLE_PARAMETERS_H

#include <RcppArmadillo.h>

using namespace Rcpp;

void res_protector(double& x);
void sample_beta_tilde(arma::mat& beta_nc_samp, arma::vec& y, arma::mat& x, arma::colvec& theta_sr, arma::vec& sig2, arma::colvec& beta_mean, int N, int d, arma::vec delta, Function Rchol); //
void sample_alpha(arma::vec& alpha_samp, arma::vec& y, arma::mat& x, arma::mat& W, arma::colvec& tau2, arma::colvec& xi2, arma::vec& sig2, arma::vec& a0, int d, Function Rchol);
void resample_alpha_diff(arma::mat& alpha_samp, arma::mat& betaenter, arma::vec& theta_sr, arma::vec& beta_mean, arma::mat& beta_diff,  arma::vec& xi2, arma::vec& tau2, int d, int N);
void sample_tau2(arma::vec& tau2_samp, arma::vec& beta_mean, double lambda2, double a_tau, int d);
void sample_xi2(arma::vec& xi2_samp, arma::vec& theta_sr, double kappa2, double a_xi, int d);
double sample_kappa2(arma::vec& xi2, double a_xi, double d1, double d2, int d);
double sample_lambda2(arma::vec& tau2, double a_tau, double e1, double e2, int d);
void sample_sigma2(arma::vec& sig2_samp, arma::vec& y, arma::mat& W, arma::vec& alpha, double c0, double C0, int N);
double sample_C0(arma::vec& sig2, double g0, double c0, double G0);
void sample_eta2(arma::vec& eta2_samp, arma::vec& alpha_samp, arma::vec& delta_samp, double d1, double d2, int d, double par_c);
void sample_p(arma::vec& delta_samp, double d1, double d2, int d, double& par_p);
void comp_prob(arma::vec& alpha_samp, arma::vec& eta2_samp, arma::vec& prob_samp, int d, double par_c, double par_p);
void sample_delta(arma::vec& prob_samp, arma::vec& delta_samp, int d);
void comp_SAVS(arma::vec& alpha_samp, arma::mat& W, arma::vec& gamma_samp, int d);
#endif

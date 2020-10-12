// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]
#include <RcppArmadillo.h>
#include <stochvol.h>
#include <progress.hpp>
#include <math.h>
#include "sample_beta_McCausland.h"
#include "sample_parameters.h"
using namespace Rcpp;

// [[Rcpp::export]]
List do_shrinkTVP_mNIG(arma::vec y,
                  arma::mat x,
                  arma::vec a0,
                  int K,
                  int niter,
                  int nburn,
                  int nthin,
                  double c0,
                  double g0,
                  double G0,
                  double eta2_d1,
                  double eta2_d2,
                  double p_d1,
                  double p_d2,
                  double par_c,
                  bool display_progress,
                  bool ret_beta_nc,
                  bool store_burn,
                  bool sv,
                  double Bsigma_sv,
                  double a0_sv,
                  double b0_sv,
                  double bmu,
                  double Bmu,
                  bool SAVS) {

  // progress bar setup
  arma::vec prog_rep_points = arma::round(arma::linspace(0, niter, 50));
  Progress p(50, display_progress);

  // Import Rs chol function
  Environment base = Environment("package:base");
  Function Rchol = base["chol"];

  // Some necessary dimensions
  int N = y.n_elem;
  int d = x.n_cols;
  int nsave;
  if (store_burn){
    nsave = std::floor(niter/nthin);
  } else {
    nsave = std::floor((niter - nburn)/nthin);
  }

  //Storage objects
  arma::cube beta_save(N+1, d, nsave, arma::fill::none);
  arma::cube sig2_save(N,1, nsave, arma::fill::none);
  arma::mat theta_sr_save(d, nsave, arma::fill::none);
  arma::mat beta_mean_save(d, nsave, arma::fill::none);
  arma::mat eta2_save(2*d, nsave, arma::fill::none);
  arma::mat prob_save(2*d, nsave, arma::fill::none);
  arma::mat delta_save(2*d, nsave, arma::fill::none);
  arma::vec p_save(nsave, arma::fill::none);
  arma::cube beta_nc_save;
  arma::mat residuals_save(N, nsave, arma::fill::none);


  // Optional Storage objects
  if (ret_beta_nc){
    beta_nc_save = arma::cube(N+1, d, nsave, arma::fill::none);
  }

  arma::vec C0_save;
  arma::vec sv_mu_save;
  arma::vec sv_phi_save;
  arma::vec sv_sigma2_save;

  if (sv == false){
    C0_save = arma::vec(nsave, arma::fill::none);
  } else {
    sv_mu_save = arma::vec(nsave, arma::fill::none);
    sv_phi_save = arma::vec(nsave, arma::fill::none);
    sv_sigma2_save = arma::vec(nsave, arma::fill::none);
  }

  // Initial values and objects
  arma::mat beta_nc_samp(d, N+1, arma::fill::none);

  arma::vec beta_mean_samp(d);
  beta_mean_samp.fill(0.1);

  arma::vec theta_sr_samp(d);
  theta_sr_samp.fill(0.2);

  arma::vec tau2_samp(d);
  tau2_samp.fill(0.1);

  arma::vec xi2_samp(d);
  xi2_samp.fill(0.1);

  arma::vec eta2_samp = arma::join_cols(xi2_samp, tau2_samp);
  arma::vec delta_samp(2*d, arma::fill::ones);
  arma::vec prob_samp(2*d, arma::fill::ones);

  double p_samp = 1;

  arma::vec h_samp(N, arma::fill::zeros);

  arma::vec alpha_samp(2*d, arma::fill::ones);

  arma::vec sig2_samp = arma::exp(h_samp);

  arma::vec residuals_samp(N, arma::fill::zeros);

  arma::vec gamma_samp(2*d, arma::fill::ones);
  arma::vec beta_mean_samp_g;
  arma::vec theta_sr_samp_g;
  arma::mat beta_nc_samp_g;
  arma::vec residuals_samp_g;

  arma::mat beta;
  if(SAVS){
    beta_mean_samp_g = arma::vec(d, arma::fill::ones);
    theta_sr_samp_g = arma::vec(d, arma::fill::ones);
    beta_nc_samp_g = arma::mat(d, N+1, arma::fill::none);
    residuals_samp_g = arma::vec(N, arma::fill::zeros);
  }


  double C0_samp = 1;

  // SV quantities
  arma::vec sv_para = {-10, 0.5, 1};
  arma::mat mixprob(10, N);
  arma::vec mixprob_vec(mixprob.begin(), mixprob.n_elem, false);
  arma::ivec r(N);
  double h0 = -10;
  double B011inv         = 1e-8;
  double B022inv         = 1e-12;
  bool Gammaprior        = true;
  double MHcontrol       = -1;
  int parameterization   = 3;
  bool centered_baseline = parameterization % 2; // 1 for C, 0 for NC baseline
  int MHsteps = 2;
  bool dontupdatemu = 0;
  double cT = N/2.0;
  double C0_sv = 1.5*Bsigma_sv;
  bool truncnormal = false;
  double priorlatent0 = -1;

  // Values for LPDS
  arma::cube m_N_save(d, 1, nsave);
  arma::cube chol_C_N_inv_save(d, d, nsave);
  arma::vec m_N_samp;
  arma::mat chol_C_N_inv_samp;

  // Values to check if the sampler failed or not
  bool successful = true;
  std::string fail;
  int fail_iter;


  // Introduce additional index post_j that is used to calculate accurate storage positions in case of thinning
  int post_j = 1;

  // Begin Gibbs loop
  for (int j = 0; j < niter; j++){
    // step a)
    // sample time varying beta.tilde parameters (NC parametrization)
    try {
      //sample_beta_tilde(beta_nc_samp, y, x, theta_sr_samp, sig2_samp, beta_mean_samp, N, d, delta, Rchol); //
     sample_beta_McCausland(beta_nc_samp, y, x, theta_sr_samp, sig2_samp, beta_mean_samp, m_N_samp, chol_C_N_inv_samp, true, N, d, Rchol);
    } catch (...){
      //beta_nc_samp.fill(nanl(""));
      if (successful == true){
        fail = "sample beta_nc";
        fail_iter = j + 1;
        successful = false;
      }
    }


    // step b)
    // sample alpha
    arma::mat x_tilde = x % (beta_nc_samp.cols(1,N)).t();
    arma::mat W = arma::join_rows(x, x_tilde);

    try {
      sample_alpha(alpha_samp, y, x, W, tau2_samp, xi2_samp, sig2_samp, a0, d, Rchol);
    } catch(...){
      alpha_samp.fill(nanl(""));
      if (successful == true){
        fail = "sample alpha";
        fail_iter = j + 1;
        successful = false;
      }
    }

    // Weave back into centered parameterization
    beta_mean_samp = alpha_samp.rows(0, d-1);
    theta_sr_samp = alpha_samp.rows(d, 2*d-1);
    arma::mat beta_nc_samp_tilde = beta_nc_samp.each_col() % theta_sr_samp;
    arma::mat betaenter = beta_nc_samp_tilde.each_col() + beta_mean_samp;

    // Difference beta outside of function (for numerical stability)
    arma::mat beta_diff_pre = arma::diff(beta_nc_samp, 1, 1);
    arma::mat beta_diff =  beta_diff_pre.each_col() % theta_sr_samp;

    // step c)
    // resample alpha
    try {
      resample_alpha_diff(alpha_samp, betaenter, theta_sr_samp, beta_mean_samp, beta_diff, xi2_samp, tau2_samp, d, N);
    } catch(...) {
      alpha_samp.fill(nanl(""));
      if (successful == true){
        fail = "resample alpha";
        fail_iter = j + 1;
        successful = false;
      }
    }



    // Calculate NC betas with new alpha
    beta_mean_samp = alpha_samp.rows(0, d-1);
    theta_sr_samp = alpha_samp.rows(d, 2*d-1);
    beta_nc_samp = betaenter.each_col() - beta_mean_samp;
    beta_nc_samp.each_col() /= theta_sr_samp;

    x_tilde = x % (beta_nc_samp.cols(1,N)).t();
    W = arma::join_rows(x, x_tilde);



    // step d)
    // sample the local shrinkage parameters eta2_j, prior probability par_p and sign parameter delta_j
    try{
      sample_eta2(eta2_samp, alpha_samp, delta_samp, eta2_d1, eta2_d2, d, par_c);
    } catch(...){
      eta2_samp.fill(nanl(""));
      if (successful == true){
        fail = "sample eta2";
        fail_iter = j + 1;
        successful = false;
      }
    }
    xi2_samp = eta2_samp.rows(0, d-1);
    tau2_samp = eta2_samp.rows(d, 2*d-1);
    try{
      comp_prob(alpha_samp, eta2_samp, prob_samp, d, par_c, p_samp);
    } catch(...){
      prob_samp.fill(nanl(""));
      if (successful == true){
        fail = "compute probability";
        fail_iter = j + 1;
        successful = false;
      }
    }

    try{
      sample_delta(prob_samp, delta_samp, d);
    } catch(...){
      delta_samp.fill(nanl(""));
      if (successful == true){
        fail = "sample delta";
        fail_iter = j + 1;
        successful = false;
      }
    }

    try{
      sample_p(delta_samp, p_d1, p_d2, d, p_samp);
    } catch(...){
      p_samp = (nanl(""));
      if (successful == true){
        fail = "sample p";
        fail_iter = j + 1;
        successful = false;
      }
    }


    // step e)
    // sample sigma2 from homoscedastic or SV case
    try {
      if (sv){
        arma::vec datastand = arma::log(arma::square(y - x * beta_mean_samp - (x % beta_nc_samp.cols(1,N).t()) * theta_sr_samp));

        arma::vec cur_h = arma::log(sig2_samp);
        stochvol::update_sv(datastand, sv_para, cur_h, h0, mixprob_vec, r, centered_baseline, C0_sv, cT,
                            Bsigma_sv, a0_sv, b0_sv, bmu, Bmu, B011inv, B022inv, Gammaprior,
                            truncnormal, MHcontrol, MHsteps, parameterization, dontupdatemu, priorlatent0);

        sig2_samp = arma::exp(cur_h);
      } else {
        sample_sigma2(sig2_samp, y, W, alpha_samp, c0, C0_samp, N);
      }
    } catch(...) {
      sig2_samp.fill(nanl(""));
      if (successful == true){
        fail = "sample sigma2";
        fail_iter = j + 1;
        successful = false;
      }
    }

    if(sv == false){
      try {
        C0_samp = sample_C0(sig2_samp, g0, c0, G0);
      } catch(...) {
        C0_samp = nanl("");
        if (successful == true){
          fail = "sample C0";
          fail_iter = j + 1;
          successful = false;
        }
      }
    }

    // adjust start of storage point depending on store_burn
    int nburn_new = nburn;
    if(store_burn){
      nburn_new = 0;
    }


    // Increment index i if burn-in period is over
    if (j > nburn_new){
      post_j++;
    }


    if (SAVS){
      comp_SAVS(alpha_samp, W, gamma_samp, d);
      beta_mean_samp_g = gamma_samp.rows(0, d-1);
      theta_sr_samp_g = gamma_samp.rows(d, 2*d-1);
      //arma::mat tmp = x % beta_nc_samp.cols(1,N).t();
      residuals_samp_g = y - x * beta_mean_samp_g - (x % beta_nc_samp.cols(1,N).t()) * theta_sr_samp_g;
      //residuals_samp_g = y - x.cols(0, K-1) * beta_mean_samp_g.rows(0, K-1) - tmp.cols(0, K-1) * theta_sr_samp_g.rows(0, K-1);
    }
    //arma::mat tmp = x % beta_nc_samp.cols(1,N).t();
    //residuals_samp = y - x.cols(0, K-1) * beta_mean_samp.rows(0, K-1) - tmp.cols(0, K-1) * theta_sr_samp.rows(0, K-1);
    residuals_samp = y - x * beta_mean_samp - (x % beta_nc_samp.cols(1,N).t()) * theta_sr_samp;
    // Store everything
    if ((post_j % nthin == 0) && (j >= nburn_new)){
      // Caluclate beta
      // This is in the if condition to save unnecessary computations if beta is not saved
      if(SAVS){
        beta =  (beta_nc_samp.each_col() % theta_sr_samp_g).each_col() + beta_mean_samp_g;
        theta_sr_save.col((post_j-1)/nthin) = theta_sr_samp_g;
        beta_mean_save.col((post_j-1)/nthin) = beta_mean_samp_g;
        residuals_save.col((post_j-1)/nthin) = residuals_samp_g;
      }else{
        beta =  (beta_nc_samp.each_col() % theta_sr_samp).each_col() + beta_mean_samp;
        theta_sr_save.col((post_j-1)/nthin) = theta_sr_samp;
        beta_mean_save.col((post_j-1)/nthin) = beta_mean_samp;
        residuals_save.col((post_j-1)/nthin) = residuals_samp;
      }


      sig2_save.slice((post_j-1)/nthin) = sig2_samp;

      beta_save.slice((post_j-1)/nthin) = beta.t();
      eta2_save.col((post_j-1)/nthin) = eta2_samp;
      p_save((post_j-1)/nthin) = p_samp;
      prob_save.col((post_j-1)/nthin) = prob_samp;
      delta_save.col((post_j-1)/nthin) = delta_samp;



      m_N_save.slice((post_j-1)/nthin) = m_N_samp;
      chol_C_N_inv_save.slice((post_j-1)/nthin) = chol_C_N_inv_samp;

      //conditional storing
      if (ret_beta_nc){
        beta_nc_save.slice((post_j-1)/nthin) = beta_nc_samp.t();
      }

      if (sv == false){
        C0_save((post_j-1)/nthin) = C0_samp;
      } else {
        sv_mu_save((post_j-1)/nthin) = sv_para(0);
        sv_phi_save((post_j-1)/nthin) = sv_para(1);
        sv_sigma2_save((post_j-1)/nthin) = sv_para(2);
      }
    }

    // Random sign switch
    for (int i = 0; i < d; i++){
      if(R::runif(0,1) > 0.5){
        theta_sr_samp(i) = -theta_sr_samp(i);
      }
    }

    // Increment progress bar
    if (arma::any(prog_rep_points == j)){
      p.increment();
    }

    // Check for user interrupts
    if (j % 500 == 0){
      Rcpp::checkUserInterrupt();
    }

    // Break loop if successful is false
    if (!successful){
      break;
    }
  }
  // return everything as a nested list (due to size restrictions on Rcpp::Lists)
  return List::create(_["sigma2"] = sig2_save,
                      _["theta_sr"] = theta_sr_save.t(),
                      _["beta_mean"] = beta_mean_save.t(),
                      _["beta_nc"] = beta_nc_save,
                      _["beta"] = beta_save,
                      _["eta2"] = eta2_save.t(),
                      _["par_p"] = p_save,
                      _["prob"] = prob_save.t(),
                      _["delta"] = delta_save.t(),
                      _["C0"] = C0_save,
                      _["residuals"] = residuals_save.t(),
                      _["sv_mu"] = sv_mu_save,
                      _["sv_phi"] = sv_phi_save,
                      _["sv_sigma2"] = sv_sigma2_save,
                      _["LPDS_comp"] = List::create(
                      _["m_N"] = m_N_save,
                      _["chol_C_N_inv"] = chol_C_N_inv_save),
                      _["success_vals"] = List::create(
                      _["success"] = successful,
                      _["fail"] = fail,
                      _["fail_iter"] = fail_iter
                      ));
}


########################################################################
#### Multivariate parcor model with shrinkage prior of inverse gamma
########################################################################
build_model <- function(F1_fwd, F1_bwd, P, m, K, n_Total){
  ### the dimension of F1_fwd and F1_bwd: n_Total * K
  ###

  ## the forward PARCOR
  n_1 <- m + 1
  n_T <- n_Total
  data <- rep(list(NA), 2*K)
  tmp_x <- F1_bwd[(n_1-m):(n_T-m), , drop = FALSE]
  y <- F1_fwd[(n_1):n_T, , drop = FALSE]
  index <- 0
  for(k in 1:K){
    index <- index + 1
    if(k==1){
      data[[index]] <- data.frame(y = y[, k, drop = FALSE], tmp_x)
    }else{
      data[[index]] <- data.frame(y = y[, k, drop = FALSE],
                                  cbind(tmp_x,
                                        y[, 1:(k-1)]))
    }
  }

  ## the backward PARCOR
  n_1 <- 1
  n_T <- n_Total - m
  tmp_x <- F1_fwd[(n_1+m):(n_T+m), , drop = FALSE]
  y <- F1_bwd[n_1:n_T, , drop = FALSE]
  for(k in 1:K){
    index <- index + 1
    if(k == 1){
      data[[index]] <- data.frame(y = y[, k, drop = FALSE],
                                  tmp_x)
    }else{
      data[[index]] <- data.frame(y = y[, k, drop = FALSE],
                                  cbind(tmp_x,
                                        y[, 1:(k-1)]))
    }
  }
  return(data)
}



PARCOR_shrinkage <- function(Y, P,
                             cpus = 1,
                             niter = 10000,
                             nburn = round(niter / 2),
                             nthin = 1,
                             hyperprior_param,
                             display_progress = TRUE,
                             ret_beta_nc = FALSE,
                             SAVS = TRUE,
                             sv = FALSE,
                             sv_param,
                             delta, S_0,
                             sample_size = 500,
                             chains = 5, DIC = FALSE, uncertainty = FALSE){
  nsave <- (niter - nburn)/nthin
  default_hyper <- list(c0 = 2.5,
                        g0 = 5,
                        G0 = 5 / (2.5 - 1),
                        eta2_d1 = 5,
                        eta2_d2 = 4,
                        p_d1 = 1,
                        p_d2 = 1,
                        par_c = 2.5/(10^5))

  if (missing(hyperprior_param)){
    hyperprior_param <- default_hyper
  }
  default_hyper_sv <- list(Bsigma_sv = 1,
                           a0_sv = 5,
                           b0_sv = 1.5,
                           bmu = 0,
                           Bmu = 1)
  if (missing(sv_param) | sv == FALSE){
    sv_param <- default_hyper_sv
  }
  ## number of time series
  K <- ncol(Y)
  ## number of time points
  n_Total <- nrow(Y)
  ### storage results
  ### stage m = 1
  result1 <- run_parcor_parallel(F1 = t(Y), delta = delta, P = 1, S_0 = S_0, sample_size = sample_size,
                                 chains = chains, DIC = DIC, uncertainty = !uncertainty)
  if(uncertainty){
    PHI_fwd_samp <- array(NA, dim = c(K^2, n_Total, P, nsave))
    PHI_bwd_samp <- array(NA, dim = c(K^2, n_Total, P, nsave))
    phi_fwd <- matrix(0, ncol = n_Total, nrow = K^2)
    phi_bwd <- matrix(0, ncol = n_Total, nrow = K^2)
    for(i in 1:nsave){
      for(j in (1+1):(n_Total-1)){
        phi_fwd[, j] <- rmvn(n = 1, mu = as.vector(result1$phi_fwd[, j, 1]),
                             sigma = result1$Cnt_fwd[[1]][[j]])
        phi_bwd[, j] <- rmvn(n = 1, mu = as.vector(result1$phi_bwd[, j, 1]),
                             sigma = result1$Cnt_bwd[[1]][[j]])
      }
      PHI_fwd_samp[, , 1, i] <- phi_fwd
      PHI_bwd_samp[, , 1, i] <- phi_bwd
    }
  }else{
    PHI_fwd <- array(NA, dim = c(K^2, n_Total, P))
    PHI_bwd <- array(NA, dim = c(K^2, n_Total, P))
    PHI_star_fwd <- array(NA, dim = c(K^2, n_Total, P))
    PHI_star_bwd <- array(NA, dim = c(K^2, n_Total, P))
    u_inv_fwd <- array(0, dim = c((K^2-K)/2, n_Total, P))
    u_inv_bwd <- array(0, dim = c((K^2-K)/2, n_Total, P))
    theta_sr <- rep(list(NA), P)
    result <- list(PHI_fwd = PHI_fwd, PHI_bwd = PHI_bwd,
                   PHI_star_fwd = PHI_star_fwd, PHI_star_bwd = PHI_star_bwd,
                   u_inv_fwd = u_inv_fwd, u_inv_bwd = u_inv_bwd,
                   SIGMA = NA, theta_sr = theta_sr)
    result$PHI_fwd[, , 1] <- result1$phi_fwd
    result$PHI_bwd[, , 1] <- result1$phi_bwd
  }


  ### stage m = 2
  data <- build_model(F1_fwd = t(result1$F1_fwd), F1_bwd = t(result1$F1_bwd), P = P,
                      m = 2, K = K, n_Total = n_Total)

  sfInit(parallel = TRUE, cpus = cpus, type = "SOCK")

  sfLibrary(PARCORwMNIG)
  sfExport("K", "niter", "nburn", "nthin", "hyperprior_param", "ret_beta_nc",
           "SAVS", "sv", "sv_param")
  for(i in 2:P){
    res_tmp <- sfLapply(1:length(data), function(x) shrinkTVP_mNIG(formula = y ~ .-1,
                                                                   data = data[[x]], K = K,
                                                                   niter = niter,
                                                                   nburn = nburn,
                                                                   nthin = nthin,
                                                                   hyperprior_param = hyperprior_param[[i]],
                                                                   display_progress = FALSE,
                                                                   ret_beta_nc = ret_beta_nc,
                                                                   SAVS = SAVS,
                                                                   sv = sv,
                                                                   sv_param = sv_param,
                                                                   simple = !uncertainty))


    if(!uncertainty){
      ### obtain median
      res <- unpack_res(res = res_tmp, m = i, n_Total = n_Total, K = K, type = "beta_median")
      data <- build_model(F1_fwd = res$resid_fwd,
                          F1_bwd = res$resid_bwd, P = P, m = i+1,
                          K = K, n_Total = n_Total)
      result$PHI_fwd[, , i] <- res$phi_fwd
      result$PHI_bwd[, , i] <- res$phi_bwd
      result$PHI_star_fwd[, , i] <- res$phi_star_fwd
      result$PHI_star_bwd[, , i] <- res$phi_star_bwd
      result$u_inv_fwd[, , i] <- res$u_inv_fwd
      result$u_inv_bwd[, , i] <- res$u_inv_bwd
      #result_all$median$theta_sr[[i]] <- res_tmp$theta_sr
      if(i == P)
        result$SIGMA <- res$SIGMA
    }else{
      res <- unpack_res_sample(res = res_tmp, m = i, n_Total = n_Total, K = K, nsave = nsave, P = P)
      data <- build_model(F1_fwd = res$resid_fwd,
                          F1_bwd = res$resid_bwd, P = P, m = i+1,
                          K = K, n_Total = n_Total)
      PHI_fwd_samp[, , i, ] <- res$phi_fwd
      PHI_bwd_samp[, , i, ] <- res$phi_bwd
    }
    cat("Stage: ", i, "/", P, "\n")
  }
  sfStop()
  if(uncertainty){
    result <- list(phi_fwd = PHI_fwd_samp, phi_bwd = PHI_bwd_samp, SIGMA = res$qSIGMA)
  }

  return(result)
}


unpack_res_sample <- function(res, m, n_Total, K, nsave, P){
  ## unpack PARCOR coefficients
  ### forward time index
  n_1_fwd <- m + 1
  n_T_fwd <- n_Total

  ### backward time index
  n_1_bwd <- 1
  n_T_bwd <- n_Total - m

  ## unpack residuals
  resid_tmp <- array(unlist(lapply(res, function(x) x$residuals)), dim = c(n_T_fwd - n_1_fwd + 1, nsave, K^2))
  resid_tmp <- aperm(resid_tmp, perm = c(1, 3, 2))

  ## unpack observational innovation
  tmp_sigma2 <- simplify2array(lapply(res, function(x) x$sigma2))
  sigma2 <- tmp_sigma2[, 1, 1:K]

  ## storage
  resid_fwd <- array(NA, dim = c(n_Total, K, nsave))
  resid_bwd <- array(NA, dim = c(n_Total, K, nsave))
  phi_fwd <- array(NA, dim = c(K^2, n_Total, nsave))
  phi_bwd <- array(NA, dim = c(K^2, n_Total, nsave))
  phi_star_fwd <- array(NA, dim = c(K^2, n_Total, nsave))
  phi_star_bwd <- array(NA, dim = c(K^2, n_Total, nsave))
  u_inv_fwd <- array(0, dim = c((K^2-K)/2, n_Total, nsave))
  u_inv_bwd <- array(0, dim = c((K^2-K)/2, n_Total, nsave))
  index_fwd <- 1
  index_bwd <- 1
  SIGMA <- array(NA, dim = c(K, K, n_T_fwd - n_1_fwd + 1, nsave))


  ### unpack the results
  for(i in 1:nsave){
    for(k in 1:K){
      ### forward
      phi_star_fwd[((k-1)*K+1):(k*K), n_1_fwd:n_T_fwd, i] <- t(res[[k]][["beta"]][-1, 1:K, i])
      if(k > 1){
        n <- nrow(t(res[[k]][["beta"]][, 1:(k-1)+K], i))
        u_inv_fwd[index_fwd:(index_fwd + n - 1), n_1_fwd:n_T_fwd, i] <- t(res[[k]][["beta"]][-1, 1:(k-1)+K, i])
        index_fwd <- index_fwd + n
      }
      ### backward
      phi_star_bwd[((k-1)*K+1):(k*K), n_1_bwd:n_T_bwd, i] <- t(res[[k+K]][["beta"]][-1, 1:K, i])

      if(k > 1){
        n <- nrow(t(res[[k+K]][["beta"]][, 1:(k-1)+K, i]))
        u_inv_bwd[index_bwd:(index_bwd + n - 1), n_1_bwd:n_T_bwd, i] <- t(res[[k+K]][["beta"]][-1, 1:(k-1)+K, i])
        index_bwd <- index_bwd + n
      }
    }


    index <- 1
    for(j in n_1_fwd:n_T_fwd){
      u_inv <- diag(K)
      u_inv[lower.tri(u_inv)] <- -u_inv_fwd[, j, i]
      u <- solve(u_inv, diag(K))
      phi_fwd[, j, i] <- as.vector(t(u%*%matrix(phi_star_fwd[, j, i], nrow = K, ncol = K, byrow = TRUE)))
      resid_fwd[j, , i] <- as.vector(u%*% matrix(resid_tmp[index, 1:K, i], ncol = 1))

      ### transformation for innovations
      if(m == P){
        SIGMA[, , j, i] <- u %*% diag(sigma2[i, ])%*%t(u)
      }
      index <- index + 1
    }
    index <- 1
    ### transformation for backward residuals
    for(j in n_1_bwd:n_T_bwd){
      u_inv <- diag(K)
      u_inv[lower.tri(u_inv)] <- -u_inv_bwd[, j, i]
      u <- solve(u_inv, diag(K))
      phi_bwd[, j, i] <- as.vector(t(u%*%matrix(phi_star_bwd[, j, i], nrow = K, ncol = K, byrow = TRUE)))
      resid_bwd[j, , i] <- as.vector(u%*% matrix(resid_tmp[index, (K+1):(2*K), i], ncol = 1))
      index <- index + 1
    }
  }
  return(list(phi_fwd = phi_fwd,
              phi_bwd = phi_bwd,
              resid_fwd = apply(resid_fwd, 1:2, median),
              resid_bwd = apply(resid_bwd, 1:2, median),
              SIGMA = SIGMA))
}


unpack_res <- function(res, m, n_Total, K, type){
  ## unpack residuals
  resid_tmp <- t(matrix(unlist(lapply(res, function(x) x$residuals)),
                        nrow = length(res), byrow = TRUE))

  ## unpack observational innovation
  ## dimension tmp_sigma2 is nsave * 1 * (2*K)
  tmp_sigma2 <- simplify2array(lapply(res, function(x) x$sigma2))
  #browser()
  sigma2 <- tmp_sigma2[, , 1:K]
  if(is.matrix(sigma2)){
    sigma2_mean <- apply(sigma2, 2, mean)
  }else if(is.array(sigma2)){
    sigma2_mean <- apply(sigma2, 2:3, mean)
  }else{
    sigma2_mean <- mean(sigma2)
  }

  ## storage
  resid_fwd <- matrix(NA, nrow = n_Total, ncol = K)
  resid_bwd <- matrix(NA, nrow = n_Total, ncol = K)
  phi_fwd <- matrix(NA, nrow = K^2, ncol = n_Total)
  phi_bwd <- matrix(NA, nrow = K^2, ncol = n_Total)
  phi_star_fwd <- matrix(NA, nrow = K^2, ncol = n_Total)
  phi_star_bwd <- matrix(NA, nrow = K^2, ncol = n_Total)
  u_inv_fwd <- matrix(0, nrow = (K^2-K)/2, ncol = n_Total)
  u_inv_bwd <- matrix(0, nrow = (K^2-K)/2, ncol = n_Total)
  index_fwd <- 1
  index_bwd <- 1

  ## unpack PARCOR coefficients
  ### forward time index
  n_1_fwd <- m + 1
  n_T_fwd <- n_Total

  ### backward time index
  n_1_bwd <- 1
  n_T_bwd <- n_Total - m

  ###
  for(k in 1:K){
    ### forward
    phi_star_fwd[((k-1)*K+1):(k*K), n_1_fwd:n_T_fwd] <- t(res[[k]][[type]][-1, 1:K])
    if(k > 1){
      n <- nrow(t(res[[k]][[type]][, 1:(k-1)+K]))
      u_inv_fwd[index_fwd:(index_fwd + n - 1), n_1_fwd:n_T_fwd] <- t(res[[k]][[type]][-1, 1:(k-1)+K])
      index_fwd <- index_fwd + n
    }
    ### backward
    phi_star_bwd[((k-1)*K+1):(k*K), n_1_bwd:n_T_bwd] <- t(res[[k+K]][[type]][-1, 1:K])

    if(k > 1){
      n <- nrow(t(res[[k+K]][[type]][, 1:(k-1)+K]))
      u_inv_bwd[index_bwd:(index_bwd + n - 1), n_1_bwd:n_T_bwd] <- t(res[[k+K]][[type]][-1, 1:(k-1)+K])
      index_bwd <- index_bwd + n
    }

  }

  SIGMA <- rep(list(NA), n_T_fwd - n_1_fwd + 1)
  ### transformation for forward residuals
  index <- 1
  for(i in n_1_fwd:n_T_fwd){
    u_inv <- diag(K)
    u_inv[lower.tri(u_inv)] <- -u_inv_fwd[, i]
    u <- solve(u_inv, diag(K))
    phi_fwd[, i] <- as.vector(t(u%*%matrix(phi_star_fwd[, i], nrow = K, ncol = K, byrow = TRUE)))
    resid_fwd[i, ] <- as.vector(u%*% matrix(resid_tmp[index, 1:K], ncol = 1))

    ### transformation for innovations
    if(length(sigma2_mean) > 1){
      SIGMA[[i]] <- u %*% diag(sigma2_mean)%*%t(u)
    }else{
      SIGMA[[i]] <- u %*% as.matrix(sigma2_mean)%*%t(u)
    }
    #browser()
    index <- index + 1
  }
  index <- 1
  ### transformation for backward residuals
  for(i in n_1_bwd:n_T_bwd){
    u_inv <- diag(K)
    u_inv[lower.tri(u_inv)] <- -u_inv_bwd[, i]
    u <- solve(u_inv, diag(K))
    phi_bwd[, i] <- as.vector(t(u%*%matrix(phi_star_bwd[, i], nrow = K, ncol = K, byrow = TRUE)))

    resid_bwd[i, ] <- as.vector(u%*% matrix(resid_tmp[index, (K+1):(2*K)], ncol = 1))
    index <- index + 1
  }

  return(list(phi_fwd = phi_fwd,
              phi_bwd = phi_bwd,
              phi_star_fwd = phi_star_fwd,
              phi_star_bwd = phi_star_bwd,
              u_inv_fwd = u_inv_fwd,
              u_inv_bwd = u_inv_bwd,
              resid_fwd = resid_fwd,
              resid_bwd = resid_bwd,
              SIGMA = SIGMA))
}




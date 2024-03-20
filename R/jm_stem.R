#' Fit joint model of longitudinal and survival data using Stochastic EM algorithm
#'
#' @import magrittr
#' @import dplyr
#' @param model1 NLME model for longitudinal data
#' @param model2 Cox PH model for survival data
#' @param covariance.model a list specifying variance-covariance matrix structure for random parameters
#' @param id a character string specifying the clusters
#' @param continuing Is it a continuous run from the previous output of jm_stem()
#' @param ini_dat Set imputed longitudinal and survival data from previous output as initial datasets
#' @param ini_par_list Set latest parameter estimates from previous output as initial estimates
#' @param sample_Sigma variance-covariance matrix for sampling random effects
#' @param iterMax  number of iterations to run
#' @param burn_in burn-in period for MH sampling
#' @param mcem_par a list of arguments for MCEM algorithm, including thin for thinning out MC chains, 
#' mc_burn for burn-in period, and S for MC sample size. By default, thin=1, mc_burn=0, and S=1 for StEM algorithm.
#' @param mc.cores number of cores used for parallel computing
#' @export
jm_stem <- function(model1, model2, covariance.model, id,
                    continuing = F, ini_dat = NULL,
                    ini_par_list = NULL, sample_Sigma=NULL, iterMax = 2500, 
                    burn_in = 50, mcem_par=list(thin = 1, mc_burn=0, S = 1), 
                    mc.cores = 1){

  memList <- list(model1, model2)
  
  # Get log likelihood of response varaibles
  RespLog <- lapply(memList, function(x){get_log_density(x)$log_density})
  random_effects <- lapply(memList, function(md){
    sapply(md$fixed$equation, function(x){strsplit(as.character(x), "~",  fixed=T)[[2]]})
  })
  fixed_params <- lapply(memList, function(md){md$fixed$start_values}) %>% unlist()
  names(fixed_params) <- lapply(memList, function(md){md$fixed$names}) %>% unlist()
  disp_params <- lapply(memList, function(md){md$sigma}) %>% unlist()
  names(disp_params) <- lapply(memList, function(md){
    if(!is.null(md$sigma)) paste0(md$response, "_sigma") else NULL }) %>% unlist()
  if(continuing==F){
    par_list <- as.list(c(fixed_params, disp_params))
  } else {
    par_list <- ini_par_list
  }

  cox_fm <- paste("survival::Surv(", model2$response, ",", model2$event, ") ~",
                  paste(model2$fixed$covariates, collapse="+")) %>% as.formula()

  if(!is.null(covariance.model$Sigma)){
    Sigma <- covariance.model$Sigma
    q <- ncol(Sigma)
  } else {
    q <- length(covariance.model$par)
    Sigma <- matrix(0, q, q)
    diag(Sigma) <- 1
    Sigma <- Sigma * covariance.model$model
  }

  fixed_names1 <- unlist(model1$fixed$names)
  p1 <- length(fixed_names1);
  p2 <- length(model2$fixed$names)
  uniqueID <- unique(c(model1$data[[id]], model2$data[[id]]))
  missing_vars <- unlist(random_effects)
  
  # set initial values for missing data
  if(continuing == F){
    missing_dat0 <- Vassign(missing_vars, rep(0, length(missing_vars)))
    dat1 <- cbind(model1$data, missing_dat0)
    dat2 <- cbind(model2$data, missing_dat0, hazard0 = 1e-16, cum_hazard0 = 1e-16)
    theta_str <- calculate_fixeff(dat1, id, model1$fixed$equation,
                                  unlist(par_list[fixed_names1]), random_effects[[1]])
    new_dat1 <- update_columns(theta_str, id, dat1)[[1]]
    new_dat2 <- update_columns(theta_str, id, dat2)[[1]]
    dup_dat1 <- new_dat1; dup_dat2 <- new_dat2
  } else {
    dup_dat1 <- ini_dat[[1]]
    dup_dat2 <- ini_dat[[2]]
  }
  
  ########################################
  ##########  Stochastic EM  #############
  ########################################
  hist_par <- list()
  hist_Ic <- list()
  hist_Sigma <- list()
  hist_basehaz <- list()
  hist_theta <- list()
  hist_loglike <- c()
  
  index1 <- split(1:nrow(dup_dat1), dup_dat1[[id]])
  index2 <- split(1:nrow(dup_dat2), dup_dat2[[id]])
  n <- length(uniqueID)
  loc <- which(covariance.model$flag==1)
  loc1 <- which(covariance.model$flag==0)
  if(is.null(sample_Sigma)){
    sample_Sigma <- diag(0.1, length(loc), length(loc))
  }

  # covariate matrix of random parameters
  X1i_list <- lapply(uniqueID, function(iid){
    lapply(model1$fixed$equation[loc], function(md){
      model.matrix(md, data= dup_dat1[dup_dat1[[id]]==iid, ][1, ])
    }) %>% Matrix::bdiag()
  })
  # covariate matrix of fixed parameters
  X2i_list <- lapply(uniqueID, function(iid){
    lapply(model1$fixed$equation[loc1], function(md){
      model.matrix(md, data= dup_dat1[dup_dat1[[id]]==iid, ][1, ])
    }) %>% Matrix::bdiag()
  })
  all_cases <- expand.grid(1:n, 1:mcem_par$S)
  
  cat("Starting (stochastic) EM algorithm...\n")
  
  m = 1
  while(m < (iterMax+1)){
    ##########################################
    ##### E-step based on Gibbs sampling #####
    ##########################################
    # (1) Sample random effects in model1 using MH method for each individual respectively
    sampling_output <- MC_sampler(posterior_theta, par_list, Sigma, sample_Sigma,
                                  dup_dat1, dup_dat2,
                                  index1, index2,
                                  fixed_names1,
                                  model1, model2,
                                  random_effects,
                                  RespLog, uniqueID,
                                  missing_vars, id, covariance.model, X1i_list,
                                  betahat = unlist(par_list[unlist(model1$fixed$names[loc])]),
                                  n, loc, mcem_par$S, mcem_par$thin, burn_in, mc.cores, 
                                  mc_burn=mcem_par$mc_burn)
    cat("sampling is done...\n")

    ##########################################
    ##### M-step for updating parameters #####
    ##########################################
    new_par <- update_parameters(sampling_output, cox_fm,
                                 missing_vars, loc,
                                 model1, model2, Sigma,
                                 mcem_par$S, n, id, uniqueID, RespLog,
                                 fixed_names1, 
                                 random_effects, covariance.model, X1i_list, 
                                 X2i_list, all_cases,
                                 par_list, mc.cores, sigma_par=names(disp_params))
    cat("par update is done...\n")
    
    # update parameters
    new_par_list <- new_par$par
    Sigma <- new_par$sigma
    info_data <- new_par$info_data

    max_diff <- mean(abs((as.numeric(par_list) - as.numeric(new_par_list))/as.numeric(par_list)))
    par_list <- new_par_list
    dup_dat1 <- info_data[[mcem_par$S]][[1]]
    dup_dat2 <- info_data[[mcem_par$S]][[2]]

    loglike = sapply(1:mcem_par$S, function(s){
      sum(with(info_data[[s]][[1]], with(new_par_list, eval(parse(text=RespLog[[1]]))))) +
        sum(with(info_data[[s]][[2]], with(new_par_list, eval(parse(text=RespLog[[2]])))))
    }) %>% Reduce("+",.)

    cat("########### iter ", m, "\n")
    cat("par=", round(unlist(par_list),2), "\n")
    cat("diff=", round(max_diff, 2), "\n")
    cat("complete loglike=", round(loglike, 2), "\n")
    cat("Sigma=\n")
    print(round(Sigma, 2))
    
    if(mcem_par$S==1){
      hist_thetam = dup_dat2[, missing_vars]
    } else {
      hist_thetam = lapply(1:mcem_par$S, function(s){
        info_data[[s]][[2]][, missing_vars]
      }) %>% Reduce("+", .)
      hist_thetam = hist_thetam/mcem_par$S
      hist_thetam[[id]] = info_data[[mcem_par$S]][[2]][[id]]
    }

    hist_par[[m]] <- unlist(par_list)
    hist_Sigma[[m]] <- as.numeric(Sigma[loc, loc])
    hist_basehaz[[m]] <- cbind(new_par$new_basehaz, iter=m)
    hist_Ic[[m]] <- as.numeric(new_par$Ic)
    hist_theta[[m]] <- hist_thetam
    hist_loglike[m] <- loglike
    m = m+1
  }

  list(hist_par = Reduce(rbind, hist_par),
       hist_Ic = Reduce(rbind, hist_Ic),
       hist_Sigma = Reduce(rbind, hist_Sigma),
       hist_basehaz = Reduce(rbind, hist_basehaz),
       hist_theta = hist_theta,
       hist_loglike = hist_loglike,
       models = list(md1 = model1,
                     md2 = model2,
                     covariance.model = covariance.model),
       last_iter_data = list(dup_dat1, dup_dat2))
}

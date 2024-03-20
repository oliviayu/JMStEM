#' Draw MC samples using Metropolis-Hastings
#'
#' @import data.table
MC_sampler <- function(posterior_theta, par_list, Sigma,
                       sample_Sigma,
                       dup_dat1, dup_dat2,
                       index1, index2,
                       fixed_names1,
                       model1, model2,
                       random_effects,
                       RespLog, uniqueID,
                       missing_vars, id, covariance.model, X1i_list, betahat,
                       n, loc, S=1, thin=1, burn_in=1, mc.cores, 
                       stopm=100, mc_burn=0, mcpackage=c("MHadaptive", "adaptMCMC", "MCMCpack")){
  
  mcpackage = match.arg(mcpackage)
  
  info_raneff <- info_data <- list()
  raneff_names <- unlist(random_effects)[loc]
  
  for(mc_s in 1:(thin*S+mc_burn)){
    ####################################
    # Sample eta and theta using MH method
    sampled_theta <- parallel::mclapply(1:n, function(i){
      sub_dat1 <- dup_dat1[index1[[i]], ]
      sub_dat2 <- dup_dat2[index2[[i]], ]
      theta_mu <- as.numeric(X1i_list[[i]] %*% matrix(betahat))
      ini_pars <- as.numeric(sub_dat1[1, missing_vars[loc]])
      prop_sigma <- sample_Sigma
      
      sub_dat1 <- dplyr::select(sub_dat1, -tidyselect::any_of(raneff_names))
      sub_dat2 <- dplyr::select(sub_dat2, -tidyselect::any_of(raneff_names))
      
      nojump <- T; rr <- 0
      while(nojump == T & rr < stopm){
        if(mcpackage=="MHadaptive"){
          res <- Metro_Hastings(posterior_theta,
                                pars = ini_pars,
                                prop_sigma =  prop_sigma,
                                par_names = raneff_names,
                                iterations = burn_in+1, burn_in=burn_in,
                                xx_names = raneff_names,
                                par_list = par_list,
                                theta_mu = theta_mu,
                                long_loglike = RespLog[[1]],
                                surv_loglike = RespLog[[2]],
                                model1=model1, loc = loc,
                                sub_dat1=sub_dat1, sub_dat2=sub_dat2,
                                Sigma=Sigma,
                                quiet =T)
          new_ranef <- tail(res$trace, 1)
        } else if (mcpackage=="adaptMCMC"){
          res <- adaptMCMC::MCMC(posterior_theta, n = burn_in + 1,
                                 init=ini_pars,
                                 scale = prop_sigma,
                                 adapt = F,
                                 xx_names = raneff_names,
                                 par_list = par_list,
                                 theta_mu = theta_mu,
                                 long_loglike = RespLog[[1]],
                                 surv_loglike = RespLog[[2]],
                                 model1=model1, loc = loc,
                                 sub_dat1=sub_dat1, sub_dat2=sub_dat2,
                                 Sigma=Sigma, showProgressBar=F)
          new_ranef <- tail(res$samples, 1)
        } else if (mcpackage == "MCMCpack"){
          res <- MCMCpack::MCMCmetrop1R(posterior_theta,
                                        theta.init=ini_pars,
                                        burnin=1,
                                        mcmc=1+burn_in,
                                        V = prop_sigma,
                                        xx_names = raneff_names,
                                        par_list = par_list,
                                        theta_mu = theta_mu,
                                        long_loglike = RespLog[[1]],
                                        surv_loglike = RespLog[[2]],
                                        model1=model1, loc = loc,
                                        sub_dat1=sub_dat1, sub_dat2=sub_dat2,
                                        Sigma=Sigma, verbose=0)
          new_ranef <- res[1,]
        }
        
        nojump <- all(new_ranef == ini_pars)
        rr <- rr+1
      }
      if(rr >= stopm){
        cat("No jump in MH algorithm...\n")
      }
      
      as.numeric(new_ranef)
    }, mc.cores = mc.cores) %>% do.call(rbind, .) %>% as.data.frame()
    
    colnames(sampled_theta) <- raneff_names
    sampled_theta[[id]] <- uniqueID
    
    # update random effects in data
    new_dat <- update_columns(sampled_theta, id, dup_dat1, dup_dat2)
    dup_dat1 <- new_dat[[1]]
    dup_dat2 <- new_dat[[2]]
    
    # (3) sample left-censored data
    censored <- (sum(dup_dat1[[model1$missing$indicator]]) > 0)
    if(censored){
      ati_ymu <- with(dup_dat1, eval(parse(text = model1$reg_equation)))
      sampled_atiy <- truncnorm::rtruncnorm(1, a=-Inf, b=model1$missing$lloq,
                                            mean = ati_ymu[dup_dat1[[model1$missing$indicator]]==1],
                                            sd = par_list[[paste0(model1$response, "_sigma")]])
      dup_dat1[[model1$response]][dup_dat1[[model1$missing$indicator]] == 1] <- sampled_atiy
    }
    
    ##########################################
    ##### info for M-step ####################
    ##########################################
    if(mc_s > mc_burn){
      info_raneff[[mc_s-mc_burn]] <- sampled_theta
      info_data[[mc_s-mc_burn]] <- list(dup_dat1, dup_dat2)
    }
  }
  
  if(thin > 1){
    info_raneff <- lapply(seq(1, thin*S, by=thin), function(i){info_raneff[[i]]})
    info_data <- lapply(seq(1, thin*S, by=thin), function(i){info_data[[i]]})
  }
  
  list(info_raneff=info_raneff,
       info_data=info_data)
}

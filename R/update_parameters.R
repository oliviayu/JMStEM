#' Update parameter estimates
#'
update_parameters <- function(sampling_output, cox_fm,
                              missing_vars, loc,
                              model1, model2, Sigma,
                              S, n, id, uniqueID, RespLog,
                              fixed_names1, 
                              random_effects, covariance.model, 
                              X1i_list, X2i_list, all_cases,
                              par_list, mc.cores, sigma_par){
  info_raneff <- sampling_output$info_raneff
  info_data <- sampling_output$info_data
  loc1 <- which(model1$random$flag==0)
  xx_names0 <- unlist(model1$fixed$names[loc])
  xx_names1 <- unlist(model1$fixed$names[loc1])
  xx_all <- names(unlist(par_list))
  p1 = length(c(xx_names0, xx_names1))
  p = p1 + length(model2$fixed$names)
  p11 = ncol(X1i_list[[1]])
  q1=length(loc)
  q11 = sum(1:length(loc))
  L2 = strMat(q1)
  
  # estimate parameters in nonlinear mixed effect models
  new_nlme_pars <- est_nlme_pars(Sigma, info_data, id, RespLog,
                                 model1, n,
                                 fixed_names1, missing_vars, uniqueID,
                                 S, random_effects, covariance.model,
                                 X1i_list, X2i_list, all_cases,
                                 par_list, mc.cores)
  new_par_list <- new_nlme_pars$beta
  fix_theta <- new_nlme_pars$fix_theta

  # estimate parameters for survival model
  dat2_all <- lapply(1:S, function(s){
    info_data[[s]][[2]]
    }) %>% do.call(rbind, .) %>% as.data.frame()
  # update fixed thetas
  dat2_all <- update_columns(fix_theta, id, dat2_all)[[1]]
  fitCOX <- survival::coxph(cox_fm, data = dat2_all, model = T)
  new_par_list[model2$fixed$names] <- as.numeric(coef(fitCOX))
  new_basehaz <- est_base_hazard(dat2_all, model2, unlist(new_par_list), Bi=dat2_all[, c(id, missing_vars)])
  new_basehaz <- new_basehaz$hazard_est[, c("hazard0", "cum_hazard0", "time")]
  names(new_basehaz)[3] <- model2$response
  
  # update fix theta and baseline haz
  for(s in 1:S){
    info_data[[s]] <- update_columns(fix_theta, id, list_input = info_data[[s]])
    info_data[[s]][[2]] <- update_columns(new_basehaz, id=model2$response, info_data[[s]][[2]])[[1]]
  }

  all_dat1 <- lapply(1:S, function(s){
    dati <- info_data[[s]][[1]]
    dati$new_id <- paste0(dati[[id]], "_s")
    dati
  }) %>% Reduce(rbind,.)
  
  Smat1 <- try(hessian_nlmepar(new_par_list, X2i_list,
                  loglike=RespLog[[1]],
                  randeff=c(random_effects[[1]][-loc], sigma_par),
                  new_dat=all_dat1, id='new_id'), silent = T)

  # Ic for parameters in random effects
  Sigma00 = new_nlme_pars$Sigma[loc, loc]
  invSigma00 <- MASS::ginv(Sigma00)

  # hessian of loglike w.r.t beta's
  Imat00 <- matrix(0, p11, p11)
  # partial derivative of beta and L's
  Imat01 <- matrix(0, nrow=p11, ncol=q11)
  # hessian of L's
  Imat11 <- matrix(0, nrow=q11, ncol=q11)
  subdat1 = group_by_at(all_dat1, c(id, "new_id")) %>% slice(1)
  Lvalt = chol(invSigma00)
  lval = Lvalt[upper.tri(Lvalt, diag=T)]
  lval = Vassign(L2$Mpar, lval)
  
  dM_list = lapply(1:q11, function(s){
    evalMat(L2$dM[,,s], lval)
  })
  d2M_list <- lapply(1:q11, function(s){
    lapply(1:q11, function(kk){
      if(kk < s){
        d2M  = NA
      } else {
        ttz <- unlist(lapply(L2$dM[,,s], function(x){ 
          Vderiv(x, L2$Mpar[kk])}))
        d2M <- evalMat(matrix(as.character(ttz), q1, q1), lval)
      }
      d2M
    })
  })
  
  for(i in 1:n){
    x = as.matrix(X1i_list[[i]])
    Imat00 = Imat00 - t(x)%*%invSigma00%*%x
    thetas_i = as.matrix(subdat1[subdat1[[id]]==uniqueID[i], random_effects[[1]][loc]])
    pred_thetai = as.numeric(x %*% matrix(unlist(new_par_list[xx_names0])))
    eps = apply(thetas_i, 2, mean) - pred_thetai
    Imat011 = c()
    Imat111 = c()
    for(s in 1:q11){
      dM = dM_list[[s]]
      Imat011 = cbind(Imat011, t(x) %*% dM %*% eps)
      Imat11k00 = rep(0, q11)
      Imat11k = sapply(s:q11, function(kk){
        dMkj = dM_list[[kk]]
        d2M = d2M_list[[s]][[kk]]
        
        0.5*matrixcalc::matrix.trace(-Sigma00 %*% dMkj %*% Sigma00 %*% dM +
                           Sigma00 %*% d2M) -
          0.5*t(eps)%*% d2M %*% eps
      })
      Imat11k00[s:q11] = Imat11k
      Imat111 = rbind(Imat111, Imat11k00)
    }
    Imat01 = Imat01 + Imat011
    Imat11 = Imat11 + Imat111
  }
  Imat11[lower.tri(Imat11)] = Imat11[upper.tri(Imat11)]
  Imat <- rbind(cbind(Imat00, Imat01),
                cbind(t(Imat01), Imat11))
  Smat2 <- -Imat
  
  if(all(class(Smat1)!='try-error')){
    Ic_est <- Matrix::bdiag(Smat1,
                     Smat2,
                     solve(fitCOX$var))
    var_order = c(xx_names1, sigma_par, xx_names0, 
                  L2$Mpar, model2$fixed$names)
    
    match1 = match(names(new_par_list), var_order)
    match2 = setdiff(1:length(var_order), match1)
    Ic_est = Ic_est[c(match1, match2), c(match1, match2)]
  } else {
    Ic_est <- NA
  }
  
  list(par=new_par_list,
       sigma=new_nlme_pars$Sigma,
       info_data=info_data,
       new_basehaz = new_basehaz,
       Ic = Ic_est)
}

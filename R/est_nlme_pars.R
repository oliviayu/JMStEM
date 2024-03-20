#' Estimate parameters in NLME model
#' 
est_nlme_pars <- function(Sigma, info_data, id, RespLog,
                          model1, n,
                          fixed_names1, missing_vars,uniqueID,
                          S, random_effects, covariance.model, 
                          X1i_list, X2i_list, all_cases,
                          par_list, mc.cores, ptol=1e-2){

  loc <- which(model1$random$flag==1)
  loc1 <- which(model1$random$flag==0)
  xx_names0 <- unlist(model1$fixed$names[loc])
  xx_names1 <- unlist(model1$fixed$names[loc1])
  dat1_all <- parallel::mclapply(1:S, function(s){
    dati <- info_data[[s]][[1]]
    dati$siter <- s
    dati
  }, mc.cores=mc.cores) %>% do.call(rbind, .)
  thetas <- parallel::mclapply(info_data, function(dat){
    subdat <- group_by_at(dat[[1]], id) %>% slice(1)
    subdat[, random_effects[[1]][loc]]
  }, mc.cores=mc.cores) %>% do.call(rbind, .)
  invA <- solve(Sigma[loc, loc])
  
  old_par <- unlist(par_list)[xx_names0]
  par_diff <- 1; kk <- 1
  
  # update parameters in random effects
  while(par_diff > ptol & kk < 50){
    # update mean parameters in random effects
    res <- parallel::mclapply(1:nrow(all_cases), function(j){
      i <- all_cases[j, 1]
      X1i <- as.matrix(X1i_list[[i]])
      res1 <- as.matrix(thetas[j, ]) %*% invA %*% X1i
      res2 <- t(X1i)%*%invA%*%X1i
      list(res1, res2)
    }, mc.cores = mc.cores)
    mat1 <- lapply(res, function(x){x[[1]]}) %>% Reduce("+", .)
    mat2 <- lapply(res, function(x){x[[2]]}) %>% Reduce("+", .)
    new_beta1 <- t(mat1 %*% solve(mat2))

    # update cov in random effects
    theta_mui <- parallel::mclapply(1:nrow(all_cases), function(j){
      i <- all_cases[j, 1]
      as.numeric(X1i_list[[i]] %*% new_beta1)}, mc.cores = mc.cores) %>%
      do.call(rbind,.)

    new_ai <- as.matrix(thetas - theta_mui)
    new_Sigma <- parallel::mclapply(1:nrow(new_ai), function(i){
      new_ai[i,] %*% t(new_ai[i,])
    }, mc.cores = mc.cores) %>% Reduce("+", .)
    new_Sigma <- new_Sigma/nrow(new_ai)
    invA <- MASS::ginv(new_Sigma)
    Sigma[loc, loc] <- new_Sigma
    par_diff <- mean(abs(new_beta1 - old_par))
    old_par <- new_beta1
    kk <- kk+1
  }
  if(kk >=100){
    stop("fixed parameters estimation in random thetas didn't converge...")
  } else {
    par_list[xx_names0] <- old_par
  }

  ##  update betas in theta's without random effects
  fix_theta_mui <- c()
  if(length(loc1) >0){
    old_par2 <- unlist(par_list)[xx_names1]
    par_diff <- 1; kk <- 1
    while(par_diff > ptol & kk < 50){
      new_beta2 <- optim(old_par2,
                         fn=fn_fixed_par, gr=score_fixed_par,
                         method = "BFGS",
                         xx_names1=xx_names1, 
                         X2i_list=X2i_list, model1=model1,
                         loglike=RespLog[[1]],
                         randeff=random_effects[[1]][loc1], par_list=par_list,
                         new_dat=dat1_all, id=id, uniqueID=uniqueID,
                         loc1=loc1)
      
      # calculate fixed theta based on beta
      fix_theta_mui <-  lapply(1:n, function(i){
        c(as.numeric(X2i_list[[i]] %*% new_beta2$par))}) %>%
        do.call(rbind, .) %>% as.data.frame()
      names(fix_theta_mui) <- random_effects[[1]][loc1]
      fix_theta_mui[[id]] <- uniqueID
      dat1_all <- update_columns(fix_theta_mui, id, dat1_all)[[1]]

      ymu <- with(dat1_all, eval(parse(text=model1$reg_equation)))
      par_list[paste0(model1$response, "_sigma")] <- sqrt(mean((dat1_all[[model1$response]]-ymu)^2))

      par_diff <- mean(abs(new_beta2$par - old_par2))
      
      old_par2 <- new_beta2$par
      kk <- kk+1
    }
    if(kk >=100){
      stop("parameters estimation in fixed thetas didn't converge...")
    } else {
      par_list[xx_names1] <- old_par2
    }
  } else {
    # update sigma
    ymu <- with(comdat, eval(parse(text=model1$reg_equation)))
    par_list[paste0(model1$response, "_sigma")] <- sqrt(mean((comdat[[model1$response]]-ymu)^2))
  }

  list(beta= par_list,
       Sigma=Sigma,
       fix_theta = fix_theta_mui)
}

#' Gradient of negative log likelihood
score_fixed_par <- function(xx, xx_names1, X2i_list, model1,
                            loglike, randeff, par_list, new_dat, id, uniqueID,
                            loc1){
  par_list[xx_names1] <- xx
  s1 <- deriv(formula(paste("~", loglike)), randeff)
  new_etas <- lapply(X2i_list, function(xi){
    as.numeric(xi %*% as.matrix(xx))
  }) %>% do.call(rbind, .) %>% as.data.frame()
  names(new_etas) <- randeff
  new_etas[[id]] <- uniqueID
  new_dat <- update_columns(new_etas, id, new_dat)[[1]]
  
  gr_eta <- with(new_dat, with(par_list, attr(eval(s1), "gradient")))

  res <- c()
  for(k in 1:length(loc1)){
    xk <- model.matrix(model1$fixed$equation[[loc1[k]]], data= new_dat)
    res <- c(res, sapply(1:ncol(xk), function(j){sum(xk[,j]*gr_eta[,k])}))
  }
  -res
}

#' Negative log likelihood
#'
fn_fixed_par <- function(xx, xx_names1, X2i_list, model1, 
                         loglike, randeff, par_list, new_dat, id, uniqueID, loc1){
  par_list[xx_names1] <- xx
  new_etas <- lapply(X2i_list, function(xi){
    as.numeric(xi %*% as.matrix(xx))
  }) %>% do.call(rbind, .) %>% as.data.frame()
  names(new_etas) <- randeff
  new_etas[[id]] <- uniqueID
  new_dat <- update_columns(new_etas, id, new_dat)[[1]]
  fn <- with(new_dat, with(par_list, eval(parse(text=eval(loglike)))))
  -sum(fn)
}
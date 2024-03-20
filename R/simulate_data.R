#' Simulate longitudinal and survival data given specified model frameworks
#' 
#' @param n sample size
#' @param truncT time when longitudinal data are truncated (in weeks)
#' @param model1 NLME model for longitudinal data
#' @param model2 Cox PH model for survival data
#' @param covariance.model a list specifying variance-covariance matrix structure for random parameters
#' @param cumhaz0 cumulative hazard function for simulation survival data
#' @param data_frequency a vector specifying how frequent (in days) longitudinal data are collected
#' @param add_var whether adding variability in data_frequency or not 
#' @export
simulate_data <- function(n = 300, truncT=24, model1, model2, covariance.model,
                          cumhaz0 = c(rep(1.5e-4, 70), rep(5e-4, 70), rep(1.5e-4, 24*7+1-140)), 
                          data_frequency= NULL, 
                          age_par=c(41, 8), sex_par=0.93, infection_type_par=0.3,
                          nadir_cd4_par = 0.3, nnrti_based_par=0.5,
                          add_var=T){
  
  memList <- list(model1, model2)
  random_effects <- lapply(memList, function(md){
    sapply(md$fixed$equation, function(x){strsplit(as.character(x), "~",  fixed=T)[[2]]})
  })
  
  ni=truncT*7+1 
  if(is.null(data_frequency)){
    data_frequency = rep(1, ni)
  }
  
  dat1 <- data.frame(PATIENT= rep(1:n, each = ni), time = rep(0:(ni-1), n)/7,
                     age = rep(rnorm(n, age_par[1], age_par[2]), each=ni),
                     sex = rep(rbinom(n, 1, sex_par), each=ni),
                     infection_type = rep(rbinom(n, 1, infection_type_par), each=ni),
                     nadir_cd4 = rep(rbinom(n, 1, nadir_cd4_par), each=ni),
                     nnrti_based = rep(rbinom(n, 1, nnrti_based_par), each=ni))
  
  dat2 <- group_by(dat1, PATIENT) %>% slice(1) %>% ungroup() %>% dplyr::select(-time)
  par_list <- c(Vassign(unlist(model1$fixed$names), model1$fixed$start_values),
                Vassign(unlist(model2$fixed$names), model2$fixed$start_values))
  
  q <- length(covariance.model$par)
  Bi <- mvtnorm::rmvnorm(n, rep(0,q), sigma = covariance.model$Sigma) %>% as.data.frame()
  names(Bi) <- covariance.model$par
  Bi$PATIENT <- 1:n
  missing_dat0 <- Vassign(unlist(random_effects), rep(0, length(unlist(random_effects))))
  dat1 <- left_join(dat1, Bi, by="PATIENT")%>% cbind(missing_dat0)
  dat2 <- left_join(dat2, Bi, by="PATIENT")%>% cbind(missing_dat0)
  
  for(i in 1:n){
    subdat <- subset(dat1, PATIENT==i)[1,]
    X1i <- lapply(model1$fixed$equation, function(md){
      model.matrix(md, data= subdat)
    }) %>% Matrix::bdiag() %>% as.matrix()
    eta_mu <- as.numeric(X1i %*% matrix(unlist(par_list[unlist(model1$fixed$names)])) + subdat[,model1$random$names])
    
    dat1[dat1$PATIENT==i, random_effects[[1]]] <- rep(eta_mu, each=ni)
    dat1[dat1$PATIENT==i, model1$response] <- with(subset(dat1, PATIENT==i),
                                                   with(par_list, eval(parse(text=model1$reg_equation)))) + 
      rnorm(ni, 0, model1$sigma)
    
    dat2[dat2$PATIENT==i, unlist(random_effects)] <- c(eta_mu)
  }
  
  sexp <- exp(with(dat2, with(par_list, eval(parse(text= model2$reg_equation)))))
  t <- (0:(ni-1))/7
  S <- sapply(1:n, function(i){
    cumhaz <- cumhaz0*sexp[i]
    Surv <- exp(-cumhaz)
    u <- runif(1)
    min(t[1-Surv >= u]) # failure time
  }) %>% unlist()
  dat2 <- mutate(dat2, event_time = pmin(S, truncT), event = as.numeric(S <= truncT))
  
  missing_vars <- c(unlist(random_effects), covariance.model$par)
  dat1 <- left_join(dat1, dat2[, c("PATIENT", "event_time")]) %>% 
    filter(time <= event_time)
  dat1 <- trim_data(dat1, data_frequency, n, add_var) %>% 
    dplyr::select(-tidyselect::any_of(missing_vars), -event_time)
  dat1[[model1$missing$indicator]] <- as.numeric(dat1[[model1$response]] < model1$missing$lloq)
  dat1[dat1[[model1$missing$indicator]]==1, model1$response] <- model1$missing$lloq
  
  dat2 <- dplyr::select(dat2, -tidyselect::any_of(missing_vars))
  list(dat1, dat2, Bi)
}


#' Trim longitudinal data based on desired data_frequency
trim_data <- function(data, data_frequency, n, add_var=F){
  res <- data.frame()
  for(i in 1:n){
    sub_dat <- dplyr::filter(data, PATIENT == i)
    if(length(data_frequency)==1){
      keep_rows <- unique(c(seq(1, nrow(sub_dat), by = data_frequency), nrow(sub_dat))) %>% sort()
    } else {
      if(add_var){
        keep_rows <- unique(c(pmin(cumsum(c(1, turb(data_frequency, 0.5))), nrow(sub_dat)), nrow(sub_dat))) %>% sort()
      } else {
        keep_rows <- unique(c(pmin(cumsum(c(1, data_frequency)), nrow(sub_dat)), nrow(sub_dat))) %>% sort()
      }
    }
    res <- rbind(res, sub_dat[keep_rows, ])
  }
  res
}

#' Add variability to data_frequency
turb <- function(data_frequency, ratio=0.5){
  bound = floor(data_frequency*ratio)
  turb = sapply(bound, function(x){
    sample(-x:x, 1)
  })
  data_frequency + turb
}

#' Estimate the baseline hazard function in the Cox model
#'
est_base_hazard <- function(surv.data, survObject, par_val, Bi){
  par_val <- as.list(par_val)
  Etime <- surv.data[[survObject$response]]
  nblock <- length(unique(Etime)) + 2
  Twindow <- c(0, sort(unique(Etime)) + 0.001, Inf)
  delta1 <- Twindow[1:(nblock - 1)] # lower window
  delta2 <- Twindow[2:nblock]  # upper window
  hazard_est <- data.frame(delta1, delta2, hazard0 = 0, cum_hazard0 = 0,
                           T_len = delta2 - delta1,
                           time = c(sort(unique(Etime)), Inf))
  linear <- parse(text= survObject$reg_equation)

  for(i in 1:nrow(hazard_est)){
    subdat <- surv.data[Etime < delta2[i] & Etime >= delta1[i], ]
    # subdat of participants who survived up to delta1[i]
    indx <- which(Etime >= delta1[i])
    subdat2 <- surv.data[indx, ]
    subBi <- Bi[indx, ]

    numer <- sum(subdat[[survObject$event]])
    exp_term <- with(subdat2, with(par_val, with(subBi, exp(eval(linear)))))
    T_len <- pmin(Etime[indx] - delta1[i], delta2[i] - delta1[i])
    denom <-  sum(exp_term*T_len)
    h0 <- numer/denom
    if(is.nan(h0)){
      hazard_est$hazard0[i] <- 1e-16
    } else {
      hazard_est$hazard0[i] <- max(1e-16, h0)
    }
    if(i ==1 ){
      hazard_est$cum_hazard0[i] <- hazard_est$hazard0[i] * (hazard_est$time[i] - hazard_est$delta1[i])
    } else{
      hazard_est$cum_hazard0[i] <- sum((hazard_est$hazard0 * hazard_est$T_len)[1:(i - 1)]) +
        hazard_est$hazard0[i]*(hazard_est$time[i] - hazard_est$delta1[i])
    }
  }

  surv.data$hazard0 <- NA
  surv.data$cum_hazard0 <- NA
  for(i in 1:nrow(surv.data)){
    kk <- which(delta1 <= Etime[i] & delta2 > Etime[i])
    surv.data$hazard0[i] <- hazard_est$hazard0[kk]
    surv.data$cum_hazard0[i] <- hazard_est$cum_hazard0[kk]
  }
  list(surv.data = surv.data, hazard_est = hazard_est)
}


#' Fit a joint model using two-step method to obtain initial estimates.
#' Note that this function is tailored for the example data and models only.
#' 
#' @param model1 NLME model for longitudinal data
#' @param model2 Cox PH model for survival data
#' @param impute_value Imputed value for left-censored longitudinal data
#' @import saemix 
#' @import dplyr 
#' @import survival
#' @export
two_step_method <- function(model1, model2, impute_value){
  
  new_dat2 <<- model1$data
  new_dat2$log_rna[new_dat2$c==1] <- impute_value
  n_par=8
  tt <- matrix(1, 8, 8)
  tt[c(2,3,5,7), ] <- 0
  tt[,c(2,3,5,7)] <- 0
  
  saemix_dat3 <- saemix::saemixData(name.data = new_dat2, header = T,
                                    sep="", na = NA, name.group = "PATIENT",
                                    name.predictors = c("time", "infection_type",
                                                        "nadir_cd4",
                                                        "nnrti_based"),
                                    name.response = "log_rna")
  saemix_md3 <- saemix::saemixModel(
    model = function(psi, id, xi){
      time <- xi[,1]
      x1 <- xi[,2]
      x2 <- xi[,3]
      x3 <- xi[,4]
      beta10 <- psi[id, 1]
      beta11 <- psi[id, 2]
      beta12 <- psi[id, 3]
      beta20 <- psi[id, 4]
      beta21 <- psi[id, 5]
      beta30 <- psi[id, 6]
      beta40 <- psi[id, 7]
      beta50 <- psi[id, 8]
      beta1 <- beta10 + beta11*x1 + beta12*x2
      beta2 <- beta20 + beta21*x3
      ypred <- (beta1*time)/(time + exp(beta2 - exp(beta30)*time))+beta40/(1+exp(beta50*time))
      return(ypred)
    },
    psi0 = matrix(model1$fixed$start_values,
                  ncol = n_par, byrow = F,
                  dimnames = list(NULL, c(paste0("beta1", 0:2),
                                          paste0("beta2", 0:1),
                                          "beta30",
                                          "beta40", "beta50"))),
    transform.par = rep(0, 8),
    fixed.estim = rep(1, 8),
    covariance.model=  tt,
    error.model="constant")
  md3 <- try(saemix::saemix(saemix_md3, data = saemix_dat3,
                            control = saemix::saemixControl(displayProgress = FALSE, print = FALSE,
                                                            save = FALSE, save.graphs =FALSE)), silent = T)
  md3_bi <- saemix::psi(md3, type="mean") %>% scale(center = F, scale = F)
  xi <- group_by(new_dat2, by=PATIENT) %>% slice(1) %>% ungroup()
  new_dat3 <- model2$data
  new_dat3$b31 <- apply(md3_bi[,1:3] * cbind(1, xi[, c("infection_type","nadir_cd4")]), 1, sum)
  new_dat3$b33 <- md3_bi[,6]
  new_dat3$b34 <- md3_bi[,8]
  
  fitCOX <- survival::coxph(survival::Surv(event_time, event) ~ b31+b33+b34+age+sex+infection_type, 
                  data = new_dat3)
  
  fix_est <- list(md1=md3@results@fixed.effects,
                  md2=coef(fitCOX))
  fix_se <- list(md1= md3@results@se.fixed,
                 md2= summary(fitCOX)$coef[,3])
  sigma_est <- list(md1 = md3@results@respar[1])
  
  list(est=fix_est,
       se = fix_se,
       sigma=sigma_est,
       models = list(md1=md3, md2=fitCOX))
}
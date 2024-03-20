#' Un-normalized log posterior of random effects theta_i
#'
posterior_theta <- function(xx, xx_names, par_list, long_loglike, surv_loglike,
                            model1, sub_dat1, sub_dat2, Sigma, loc, theta_mu){
  # update parameters
  sub_dat1 <- Vassign(xx_names, xx) %>%
    cbind(dplyr::select(sub_dat1, -tidyselect::any_of(xx_names)), .)
  sub_dat2 <- Vassign(xx_names, xx) %>%
    cbind(dplyr::select(sub_dat2, -tidyselect::any_of(xx_names)), .)

  # log likelihood of longitudinal data y|theta
  fn1 <- with(par_list, with(sub_dat1, eval(parse(text = long_loglike)))) %>% sum()
  # log likelihood of theta
  theta_res <- matrix(xx - theta_mu)
  fn2 <- -0.5*t(theta_res)%*%solve(Sigma[loc, loc])%*%theta_res
  # log likelihood of survival data, given psi(theta)
  fn3 <- with(par_list, with(sub_dat2, eval(parse(text = surv_loglike))))

  fn1+fn2+fn3
}


#' Un-normalized log posterior of random effects theta_i
#'
posterior_theta_nlme <- function(xx, xx_names, par_list, long_loglike, 
                            model1, sub_dat1, Sigma, loc, theta_mu){
  # update parameters
  sub_dat1 <- Vassign(xx_names, xx) %>%
    cbind(dplyr::select(sub_dat1, -tidyselect::any_of(xx_names)), .)

  # log likelihood of longitudinal data y|theta
  fn1 <- with(par_list, with(sub_dat1, eval(parse(text = long_loglike)))) %>% sum()
  # log likelihood of theta
  theta_res <- matrix(xx - theta_mu)
  fn2 <- -0.5*t(theta_res)%*%solve(Sigma[loc, loc])%*%theta_res

  fn1+fn2
}

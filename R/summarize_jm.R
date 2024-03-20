#' Summarize output of jm_stem()
#' 
#' @export
summarize_jm = function(resi, burn_in=500, id='id'){
  
  kk = (burn_in +1):nrow(resi$hist_par)
  mm <- length(kk)
  
  p1 = length(jm_res$models$md1$fixed$start_values) + length(jm_res$models$md2$fixed$start_values)
  mu <- apply(resi$hist_par[kk, 1:p1], 2, mean)
  ni <- length(unique(resi$models$md1$data[[id]]))
  
  if(length(resi$models)>2){
    aaa <- t(matrix(resi$hist_basehaz$hazard0, nrow = 
                      length(unique(resi$hist_basehaz$event_time))))
    aaa <- aaa[, 1:(ncol(aaa)-2)]
  } else {
    aaa <- NULL
  }
  
  q = sqrt(length(resi$hist_Sigma[1,]))
  p2 = length(jm_res$hist_par[1,]) + sum(1:q) 
  ## complete data info: within variance
  Ic = apply(resi$hist_Ic[kk, ], 2, mean)/ni %>% matrix(p2, p2)
  Wvar <- solve(Ic)
  Lmat = lapply(kk, function(i){
    Lmat = chol(solve(matrix(resi$hist_Sigma[i,], q, q)))
    as.numeric(Lmat[upper.tri(Lmat, diag=T)])
  }) %>% do.call(rbind, .)
  all_est = cbind(resi$hist_par[kk, ], Lmat)
  ### sample covariance matrix
  Bvar00 <- cov(all_est)*ni
  Bvar <- cov(all_est[seq(1, nrow(all_est), by=1), ])*ni
  
  ttt <- MASS::ginv(Bvar)%*%Wvar
  Fmat <- expm::sqrtm(ttt%*%ttt/4 + diag(1, p2, p2)) - ttt/2

  Iobs_inv = Wvar%*%MASS::ginv(diag(1, p2, p2) - Fmat)

  II <- diag(1, nrow(Fmat), nrow(Fmat))
  FF <- Fmat
  FFmm <- Fmat
  for(jj in 1:(mm-1)){
    FFmm <- FFmm%*%Fmat
  }
  tunef <- II + (II-solve(II + FF))/mm + 2/mm*(II- solve(II + FF))%*%FF%*%solve(II-FF) -
    2/mm^2*(II-solve(II+FF))%*%FF%*%(II-FFmm)%*%(solve(II-FF)%*%solve(II-FF))
  sd1 = sqrt(diag(Iobs_inv%*%tunef)/ni)[1:p1]

  threval <- 1.96
  data.frame(est= mu, sd= sd1, 
             ci_low = mu - threval*sd1,
             ci_up = mu + threval*sd1,
             p = (1-pnorm(abs(mu/sd1)))*2)
}

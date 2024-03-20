
#' @importFrom Matrix bdiag
hessian_nlmepar <- function(par_list, X2i_list, loglike, randeff, new_dat, id){
  uniqueID <- unique(new_dat[[id]])
  s1 <- deriv(formula(paste("~", loglike)), randeff, hessian=T)
  he_eta <- with(new_dat, with(par_list, attr(eval(s1), "hessian")))
  Hval <- lapply(1:length(uniqueID), function(i){
    loci <- which(new_dat[[id]] == uniqueID[i])
    Xi <- Matrix::bdiag(X2i_list[[i]], 1) %>% as.matrix()
    lapply(1:length(loci), function(j){
      t(Xi)%*%he_eta[loci[j], randeff, randeff]%*%Xi
    }) %>% Reduce("+", .)
  }) %>% Reduce("+", .)
  -as.matrix(Hval)
}


score_nlmepar <- function(par_list, model1, n,
                          loglike, randeff, new_dat, id, byall=T){
  s1 <- deriv(formula(paste("~", loglike)), randeff, hessian=F)
  gr_eta <- with(new_dat, with(par_list, attr(eval(s1), "gradient")))

  if(!is.null(model1$fixed$equation)){
    score_values <- lapply(1:length(randeff), function(k){
      xk <- model.matrix(model1$fixed$equation[[k]], data= new_dat)
      lapply(1:ncol(xk), function(j){xk[,j]*gr_eta[,k]}) %>% do.call(cbind,.)
    }) %>% do.call(cbind,.) %>% as.data.frame()
  } else {
    score_values = gr_eta
  }
  
  if(byall==T){
    res = apply(score_values, 2, sum)
  } else {
    score_values[[id]] <- new_dat[[id]]
    score_id <- group_by_at(score_values, id) %>%
      summarise(across(everything(), sum))
    res = score_id
  }
  res
}
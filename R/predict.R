#' Predict longitudinal data based on Monolix output
#' 
predict_saem <- function(theta_dat, jm_res, id = 'PATIENT', time='week', smooth=F){
  
  md1 <- jm_res$models[[1]]
  if(smooth!=T){
    dat1 <- left_join(md1$data, theta_dat, by=id)
  } else {
    dat1 <- dplyr::group_by(md1$data,across(all_of(id))) %>% slice(n()) %>%
      left_join(theta_dat, by=id)
    
    dat1 <- dat1[rep(seq_len(nrow(dat1)), dat1[[time]]), ]
    dat1[[time]] <- do.call(c, sapply(table(dat1[[id]]), function(ss){seq_len(ss)}))
  }
  
  dat1$yihat <- with(dat1, eval(parse(text=md1$reg_equation)))
  dat1 <- left_join(dplyr::select(dat1, -log_rna), 
                    md1$data[, c(id, time, "log_rna")],
                    by=c(id, time))

  dat1
}

#' Predict longitudinal data based on JMStEM output
#' 
#' @export
predict_jm <- function(jm_res, kk, id='PATIENT',
                       par_est=NULL,
                       time='time', smooth=F, method='mean'){
  md1 <- jm_res$models[[1]]
  
  if(method=='mean'){
    theta_mat <- lapply( kk, function(k){
      jm_res$hist_theta[[k]][, c(paste0("theta", 1:5))] %>%
        as.matrix()
    }) %>% Reduce("+", .)
    theta_dat <- as.data.frame(theta_mat/length(kk))
  } else if(method=='median'){
    theta_mat=lapply(paste0("theta", 1:5), function(thetai){
      mati=lapply( kk, function(k){
        jm_res$hist_theta[[k]][, thetai] %>%
          as.matrix()
      }) %>% do.call(cbind, .)
      apply(mati, 1, median)
    }) %>% do.call(cbind,.)
    theta_dat <- as.data.frame(theta_mat)
  } else if(method=='marginal'){
    if(is.null(par_est)){
      par0 = summarize_jm(jm_res, kk = kk, p=15, id=id)
      par1 = par0[,1]
      names(par1) = rownames(par0)
    } else {
      par1 = par_est
    }
    subdat= group_by(md1$data, across(all_of(id))) %>% slice(1)
    subdat[, paste0('theta',1:5)] <- 0

    theta_dat <- lapply(1:nrow(subdat), function(i){
      Xi <- lapply(md1$fixed$equation, function(md){
        model.matrix(md, data= subdat[i, ])
      }) %>% Matrix::bdiag() %>% as.matrix()
      as.numeric(Xi %*% par1[unlist(md1$fixed$names)])
    }) %>% do.call(rbind, .) %>% as.data.frame()
  }

  names(theta_dat) <- c(paste0("theta", 1:5))
  theta_dat[[id]] <- jm_res$models$md2$data[[id]]
  
  #### Re-estimate random effects using empirical Bayes method
  
  if(smooth!=T){
    dat1 <- left_join(md1$data, theta_dat, by=id)
  } else {
    dat1 <- group_by(md1$data, across(all_of(id))) %>% slice(n()) %>%
      left_join(theta_dat, by=id)
    
    dat1 <- dat1[rep(seq_len(nrow(dat1)), floor(dat1[[time]])), ]
    dat1[[time]] <- do.call(c, sapply(table(dat1[[id]]), function(ss){seq_len(ss)}))
  }
  
  dat1$yihat <- with(dat1, eval(parse(text=md1$reg_equation)))
  dat1 <- left_join(dplyr::select(dat1, -log_rna), 
                    md1$data[, c(id, time, "log_rna")],
                    by=c(id, time))
  dat1
}

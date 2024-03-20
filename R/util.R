#' Assign names to a vector
#'
Vassign <- function(name, value){
  dat <- data.frame(as.list(value))
  names(dat) <- name
  return(dat)
}


#' Evaluate a matrix given texts
#'
#' @param mat a matrix of texts
#' @param data data used for evaluation
#' @param par_val parameter values needed for evaluation
#' @param byrow logical. If TRUE (the default) the matrix is filled by rows.
evalMat <- function(mat, data = NULL, par_val = NULL, byrow = TRUE){
  q <- sqrt(length(mat))
  if(q%%1 != 0){
    stop("mat must be a square matrix.")
  }
  par_val <- as.list(par_val)
  
  res <- lapply(mat, function(x){
    if(class(x) == "character"){
      res <- eval(parse(text = x), c(data, par_val)) %>% as.numeric()
    } else {
      res <- eval(x, c(data, par_val)) %>% as.numeric()
    }
    if(!is.null(data)){
      if(nrow(data) > 1 & length(res) == 1){
        # result is independent from data values
        res <- res*nrow(data)
      }
    }
    sum(res)
  })
  
  matrix(as.numeric(res), q, q, byrow = byrow)
}


#' This function returns a matrix with elements in string.
#' Moreover, it returns L(l)'L(l) = SIGMA, which is a spherical
#' parameterization with diagonal elements equal to 1.
#'
#' @importFrom Deriv Simplify
strMat <- function(q2, method = "cholesky"){
  L <- c()
  Mpar <- c()
  
  if(method == "spherical"){
    L <- rbind(L, c(1, rep(0, q2 - 1)))
    for(i in 2:q2){
      l0 <- 1
      Li <- c()
      for(j in 1:i){
        m0 <- paste('L', i, 2:min((j+1), i), sep = "")
        Mpar <- c(Mpar, m0)
        
        if(j < i){
          sin0 <- rep(c("sin(", "cos("), c(length(m0) - 1, 1))
        } else {
          sin0 <- rep("sin(",  length(m0))
        }
        l1 <- paste(sin0, m0, rep(")", length(m0)), sep = "")
        Li <- c(Li, paste(c(l0, l1), collapse = "*"))
      }
      L <- rbind(L, c(Li, rep(0, q2 - i)))
    }
  } else if(method == "cholesky"){
    for(i in 1:q2){
      m0 <- paste0("L", i, 1:i)
      L <- rbind(L, c(m0, rep(0, q2-i)))
      Mpar <- c(Mpar, m0)
    }
  }
  
  L <- t(L)
  
  M <- matrix(NA, q2, q2)  # M=L'L
  for(i in 1:q2){
    for(j in 1:q2){
      M[i, j] <- Deriv::Simplify(paste(L[, i], L[, j], sep = "*", collapse = "+"))
    }
  }
  
  if(method == "spherical"){ diag(M) <- "1" }
  
  Mpar <- unique(Mpar)
  Mp <- length(Mpar)
  dM <- array(NA, dim = c(q2, q2, Mp))
  for(i in 1:Mp){
    ttz <- unlist(lapply(M, function(x){ Vderiv(x, Mpar[i]) }))
    dM[,,i] <- matrix(as.character(ttz), q2, q2)
  }
  
  dL <- array(NA, dim = c(q2, q2, Mp))
  for(i in 1:Mp){
    ttz <- unlist(lapply(t(L), function(x){ Vderiv(x, Mpar[i]) }))
    dL[,,i] <- matrix(as.character(ttz), q2, q2)
  }
  
  return(list(M=M, Mpar=Mpar, dM=dM, L=t(L), dL=dL))
}


#' Compute derivatives of an expression w.r.t
#' a vector of parameters respectively
#'
Vderiv <- function(lik1, pars){
  lik <- parse(text = lik1)
  q <- length(pars)
  result <- as.list(rep(NA, q))
  
  for(i in 1:q){
    result[[i]] <- Deriv::Simplify( D(lik, pars[i]) )
    names(result)[i] <- pars[i]
  }
  return(result)
}
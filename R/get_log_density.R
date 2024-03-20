#' Get log density function of (nonlinear) (generalized) linear model or survial model
#'
get_log_density <- function(model_object){
  if(model_object$modeltype %in% c("nlme", 'lme','glme')){
    get_log_density_glm(model_object)
  } else if (model_object$modeltype == "survival"){
    get_log_density_ph(model_object)
  }
}

#' Get log density function of (generalized) linear model
#'
get_log_density_glm <- function(model_object){
  if(is.null(model_object)){
    return(NULL)
  }

  additional_param <- NULL
  left_censored <- !is.null(model_object$left_censoring)

  if(!left_censored){

    if(model_object$distribution == "binomial"){

      res <- paste(model_object$response, "*(",
                   model_object$reg_equation, ")-", "log(1+exp(",
                   model_object$reg_equation, "))")

    } else if(model_object$distribution == "normal"){

      sigma <- paste0(model_object$response, "_sigma")
      res <- paste("- 0.5*(", model_object$response, "- (",
                   model_object$reg_equation, "))^2/",
                   sigma, "^2-log(", sigma, ")-0.5*log(2*pi)")
      additional_param <- c(additional_param, sigma)

    } else if(model_object$distribution == "poisson"){

      res <- paste(model_object$response, "*(",
                   model_object$reg_equation, ")-exp(",
                   model_object$reg_equation, ")-log(factorial(",
                   model_object$response, "))")
    }
  } else {

    sigma <- paste0(model_object$response, "_sigma")
    log_density <- paste("- 0.5*(", model_object$response,
                         "- (", model_object$reg_equation, "))^2/",
                         sigma, "^2-log(", sigma, ")-0.5*log(2*pi)")
    Cresp <- model_object$left_censoring$indicator

    if(model_object$left_censoring$method == "tobit"){
      # limit_val <- model_object$response
      limit_val <- model_object$left_censoring$limit_value
    } else if(model_object$left_censoring$method == "truncated"){
      stopifnot(!is.null(model_object$left_censoring$limit_value))
      limit_val <- model_object$left_censoring$limit_value
    }

    CstdNmu <- paste("(", limit_val, "-(",
                     model_object$reg_equation, "))/", sigma)
    leftint <- paste("1/(1+exp(-1.702*", CstdNmu, "))")

    if(model_object$left_censoring$method == "tobit"){
      # A Tobit model, assuming left-censored data follow normal distribution
      res <- paste( "(", log_density, ")*(1-", Cresp, ") +log(", leftint, ")*",Cresp )

    } else if(model_object$left_censoring$method == "truncated"){
      # Model uncensored data only, assumed following truncted normal distribution
      res <- paste("(", log_density, "-log(1-", leftint, "))*(1-", Cresp, ")")
    }

    additional_param <- c(additional_param, sigma)
  }

  list(log_density = paste("(", res, ")"),
       additional_param = additional_param)
}

#' Get log density function of propotional hazard model
#'
get_log_density_ph <- function(model_object){
  if(is.null(model_object$distribution)){
    # Cox PH model
    log_density <- paste(model_object$event, "*( log(hazard0) +",
                         model_object$reg_equation, ") - exp(",
                         model_object$reg_equation, ")*cum_hazard0")

    additional_param <- c("hazard0", "cum_hazard0")

  } else if(model_object$distribution == "weibull"){
    # Weibull PH model
    log_density <- paste(model_object$event, "*( Wlogscale+log(Wshape)+log(",
                         model_object$response, ")*(Wshape-1) +",
                         model_object$reg_equation, ") - exp(Wlogscale+",
                         model_object$reg_equation, ")*",
                         model_object$response, "^Wshape")

    additional_param <- c("Wlogscale", "Wshape")
  }

  list(log_density = log_density,
       additional_param = additional_param)
}


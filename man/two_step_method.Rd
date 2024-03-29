% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/two_step_method.R
\name{two_step_method}
\alias{two_step_method}
\title{Fit a joint model using two-step method to obtain initial estimates.
Note that this function is tailored for the example data and models only.}
\usage{
two_step_method(model1, model2, impute_value)
}
\arguments{
\item{model1}{NLME model for longitudinal data}

\item{model2}{Cox PH model for survival data}

\item{impute_value}{Imputed value for left-censored longitudinal data}
}
\description{
Fit a joint model using two-step method to obtain initial estimates.
Note that this function is tailored for the example data and models only.
}

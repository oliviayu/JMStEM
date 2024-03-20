library(devtools)
load_all(".")

#####################################################
################ Set-up joint models ################
#####################################################
# lower limit of quantification for longitudinal data
lloq <- log10(50) 

# An NLME model for longitudinal data with left-censoring
model1 <- list(
  modeltype = "nlme",
  response = "log_rna",
  reg_equation = "(theta1*time)/ (time + exp(theta2 - exp(theta3)*time)) + 
  theta4/(1+exp(theta5*time))",
  distribution = 'normal',
  fixed = list(equation = list(theta1 ~ 1 + infection_type + nadir_cd4,
                               theta2 ~ 1 + nnrti_based,
                               theta3 ~ 1,
                               theta4 ~ 1,
                               theta5 ~ 1),
               names = list(paste0("beta1", 0:2),
                            paste0("beta2", 0:1),
                            "beta3", "beta4", "beta5"),
               start_values = c(3.8,  -0.45, -0.35, 
                                5.5, 3, 0.75, 2, 0.2)),
  random = list(flag=c(1,1,1,0,1), names=paste0("a", 1:5)),
  sigma = 0.5,
  missing = list(indicator = "c",
                 lloq = lloq))

# A Cox PH model for survival data
model2 <- list(modeltype = "survival",
               response = "event_time",
               event = "event",
               reg_equation = "lambda1 * theta1 + 
               lambda2 * theta3 + lambda3 * theta5 + 
               lambda4*age + lambda5*sex + lambda6*infection_type",
               distribution = NULL,
               fixed = list(names = paste0("lambda",1:6),
                            start_values = c(2.5, 1.5, -9.5, 
                                             -0.03, -1, 2.3),
                            covariates = c("theta1", "theta3","theta5",
                                           "age", "sex", "infection_type")))

# Specify variance-covariance matrix in multivariate normal distribution 
# for all the random parameters theta1-theta5.
cov_flag <- matrix(1, nrow=5, ncol = 5)
cov_flag[4,] <- 0 # no randomness in theta4
cov_flag[,4] <- 0
covariance.model <- list(par = c(paste0("a", 1:5)),
                         flag = model1$random$flag,
                         model = cov_flag,
                         Sigma = NULL)
covariance.model$Sigma <- matrix(c(0.7, 0.5, 0.02, 0, 0.1,
                                   0.5, 4.5, -0.2, 0, 0.3,
                                   0.02, -0.2, 0.07, 0,  0,
                                   0, 0, 0, 0, 0,
                                   0.1, 0.3, 0, 0, 0.05), 5, 5)

###########################################################
######### Simulate longitudinal and survival data #########
###########################################################
### Example 1: constant baseline hazard with one-month delay
# haz0 = c(rep(0, 30), rep(5e-6, 139))

### Example 2: Weibull-distributed baseline hazard
truncT = 24
t = (0:(truncT*7))/7
shape=3
lambda = 1e-8
haz0 = lambda*shape*t^(shape-1)
cumhaz0 = cumsum(haz0)

set.seed(101)
sample_dat <- simulate_data(n = 300, truncT, model1, model2, 
                            covariance.model,
                            cumhaz0 = cumhaz0,
                            data_frequency = c(rep(3, 10), rep(7, 21), rep(30, 10)),
                            add_var = T)
rm_id <- names(which(table(sample_dat[[1]]$PATIENT) <= 3)); length(rm_id)
dat1 <- subset(sample_dat[[1]], !PATIENT %in% rm_id)
dat2 <- subset(sample_dat[[2]], !PATIENT %in% rm_id)
model10 <- model1
model20 <- model2
covariance.model0 <- covariance.model
model10$data <- dat1
model20$data <- dat2

# Run joint model based on two-step method to obtain initial estimates
nv_res <- two_step_method(model10, model20, impute_value = lloq/2)

# Update initial estimates
model10$fixed$start_values <- as.numeric(nv_res$est$md1)
model20$fixed$start_values <- as.numeric(nv_res$est$md2)
model10$sigma <- as.numeric(nv_res$sigma$md1)
covariance.model0$Sigma[c(1,2,3,5), c(1,2,3,5)] <- nv_res$models$md1@results@omega[c(1,4,6,8), c(1,4,6,8)]

# Fit joint model based on StEM algorithm
jm_res <- jm_stem(model10, model20, covariance.model0, id='PATIENT', iterMax = 2500)
# Summarize results
summarize_jm(jm_res, id='PATIENT', burn_in = 500)



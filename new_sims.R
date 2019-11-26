library(MASS)
library(propagate)
library(Renvlp)

X.simulator <- function(n, p, lxvar, uxvar, lcor = 0, ucor = 0, 
                        mean_zero = TRUE, lmu = NULL, umu = NULL, 
                        standardize_n = FALSE) {
  
  # Constructing means and covariance matrix
  if(mean_zero == TRUE){
    mu <- rep(0, p)
  } else {
    mu <- runif(p, lmu, umu)
  }
  
  cor_vec <- runif(p, sqrt(lcor), sqrt(ucor))
  cor_mat <- cor_vec %*% t(cor_vec)
  x_variances <- runif(p, lxvar, uxvar)
  diag(cor_mat) <- 1
  sigma <- cor2cov(cor_mat, x_variances)
  
  # Creating X matrix
  X <- mvrnorm(n, mu, sigma, empirical = T)
  if(standardize_n == FALSE){
    return(X)
  } else {
    return(scale(X))
  }
}

X.nonnorm.simulator <- function(n, bgl, bgu, logis, df, standardize_nn = F) {

  # Creating non-normal predictors
  g <- rgamma(n, shape = bgl, rate = bgu)
  b <- rbeta(n, shape1 = bgl, shape2 = bgu)
  l <- rlogis(n, scale = logis)
  t <- rt(n, df)
  
  # Creating X matrix
  X <- cbind(g, b, l, t)
  if(standardize_nn == FALSE){
    return(X)
  } else {
    return(scale(X))
  }
}

Y.simulator <- function(X, n, sd_in, beta_vec){
  # Creating Responses
  noise <- rnorm(n, mean = 0, sd = sd_in) # increased noise
  Y <- (X %*% unlist(beta_vec)) + noise
  return(Y)
}

non.norm.env.estimator <- function(X, Y, beta_vec) {
  
  # Fitting Envelopes
  dim <- u.xenv(X, Y)
  env_aic <- xenv(X, Y, u = dim$u.aic)
  env_bic <- xenv(X, Y, u = dim$u.bic)
  env_lrt <- xenv(X, Y, u = dim$u.lrt)
  
  # Computing Differences
  env_aic_diff <- env_aic$beta - beta_vec
  env_bic_diff <- env_bic$beta - beta_vec
  env_lrt_diff <- env_lrt$beta - beta_vec
  
  output <- list("aic" = env_aic_diff,
                 "bic" = env_bic_diff,
                 "lrt" = env_lrt_diff)
  
  return(output)
}


non.norm.ols.estimator <- function(X, Y, beta_vec){
  # Fitting OLS
  ols_fit <- lm(Y ~ X)
  ols_diff <- coef(ols_fit)[-1] - beta_vec
  return(ols_diff)
}

# Executes simulations with normally distributed data
norm.simulator <- function(trials, n, p, lxvar, uxvar, lcor = 0, ucor = 0, 
                           mean_zero = TRUE, lmu = NULL, umu = NULL, 
                           standardize_n = FALSE, lbeta = -2, ubeta = 2, 
                           sd_in = 15) {
  
  # Making output dataframes
  aic_output_df <- data.frame(matrix(nrow = trials, ncol = p))
  bic_output_df <- data.frame(matrix(nrow = trials, ncol = p))
  lrt_output_df <- data.frame(matrix(nrow = trials, ncol = p))
  ols_output_df <- data.frame(matrix(nrow = trials, ncol = p))
  
  for(i in 1:trials){
    # Generating Data
    X <- X.simulator(n, p, lxvar, uxvar, lcor = 0, ucor = 0, 
                           mean_zero = TRUE, lmu = NULL, umu = NULL, 
                           standardize_n = FALSE)
    beta_vec <- runif(ncol(X), lbeta, ubeta)
    Y <- Y.simulator(X, n, sd_in, beta_vec)
    
    # Storing output
    aic_output_df[i, ] <- non.norm.env.estimator(X, Y, beta_vec)$aic
    bic_output_df[i, ] <- non.norm.env.estimator(X, Y, beta_vec)$bic
    lrt_output_df[i, ] <- non.norm.env.estimator(X, Y, beta_vec)$lrt
    ols_output_df[i, ] <- non.norm.ols.estimator(X, Y, beta_vec)
    
  }
  return(list("aic" = aic_output_df,
              "bic" = bic_output_df,
              "lrt" = lrt_output_df,
              "ols" = ols_output_df))
}

# Executes simulation with non-normal distritbutions
non.norm.simulator <- function(trials, n, p, lxvar, uxvar, lcor = 0, ucor = 0, 
                               mean_zero = TRUE, lmu = NULL, umu = NULL, 
                               standardize_n = FALSE, standardize_nn = FALSE,
                               lbeta = -2, ubeta = 2, sd_in = 15, bgl, bgu, logis, df) {
  
  # Making output dataframes
  aic_output_df <- data.frame(matrix(nrow = trials, ncol = (4 + p)))
  bic_output_df <- data.frame(matrix(nrow = trials, ncol = (4 + p)))
  lrt_output_df <- data.frame(matrix(nrow = trials, ncol = (4 + p)))
  ols_output_df <- data.frame(matrix(nrow = trials, ncol = (4 + p)))
  
  for(i in 1:trials){
    # Generating Data
    X <- cbind(X.simulator(n, p, lxvar, uxvar, lcor = 0, ucor = 0, 
                           mean_zero = TRUE, lmu = NULL, umu = NULL, 
                           standardize_n = FALSE),
               X.nonnorm.simulator(n, bgl, bgu, logis, df, standardize_nn = FALSE))
    beta_vec <- runif(ncol(X), lbeta, ubeta)
    Y <- Y.simulator(X, n, sd_in, beta_vec)
    
    # Storing output
    aic_output_df[i, ] <- non.norm.env.estimator(X, Y, beta_vec)$aic
    bic_output_df[i, ] <- non.norm.env.estimator(X, Y, beta_vec)$bic
    lrt_output_df[i, ] <- non.norm.env.estimator(X, Y, beta_vec)$lrt
    ols_output_df[i, ] <- non.norm.ols.estimator(X, Y, beta_vec)
  }
  
  list_output <- list("aic" = aic_output_df,
                      "bic" = bic_output_df,
                      "lrt" = lrt_output_df,
                      "ols" = ols_output_df)
  
  return(list_output)
}

# Makes normal control lists
norm.control.lister <- function(trials, n, p, lxvar, uxvar, lcor, ucor, 
                                mean_zero = TRUE, lmu = NULL, umu = NULL, 
                                standardize_n = FALSE, lbeta = -2, 
                                ubeta = 2, sd_in){
  
  return(list("trials" = trials, "n" = n, "p" = p, "lxvar" = lxvar, 
              "uxvar" = uxvar, "lcor" = lcor, "ucor" = ucor, 
              "mean_zero" = mean_zero, "lmu" = lmu, "umu" = umu,
              "standardize_n" = standardize_n, "lbeta" = lbeta, "ubeta" = ubeta, 
              "sd_in" = sd_in))
}

# Creates non-normal control lists
non.norm.control.lister <- function(trials, n, p, lxvar, uxvar, lcor, ucor, 
                                mean_zero = TRUE, lmu = NULL, umu = NULL, 
                                standardize_n = FALSE, lbeta = -2, 
                                ubeta = 2, sd_in, bgl, bgu, logis, df){
  
  return(list("trials" = trials, "n" = n, "p" = p, "lxvar" = lxvar, 
              "uxvar" = uxvar, "lcor" = lcor, "ucor" = ucor, 
              "mean_zero" = mean_zero, "lmu" = lmu, "umu" = umu,
              "standardize_n" = standardize_n, "lbeta" = lbeta, "ubeta" = ubeta, 
              "sd_in" = sd_in, "bgl" = bgl, "bgu" = bgu, "logis" = logis, "df" = df))
}

# Normal Simulations----
## Control Vectors
nsim <- 5000
sample_sizes <- c(25, 75, 250)
parameters <- c(0.1, 0.3)
lower_vars <- c(1/10, 4)
upper_vars <- c(4, 15)
lower_cors <- c(0, 0.3, 0.7)
upper_cors<- c(0.3, 0.7, 1)
errors <- c(10, 20)

simulations <- length(sample_sizes) * length(parameters) * length(lower_cors) * length(errors) * length(lower_vars)
index <- 1
norm_control_list <- list()

## Creating a list of control lists
for(e in 1:length(errors)){
  for(c in 1:length(lower_cors)){
    for(v in 1:length(lower_vars)){
      for(p in 1:length(parameters)){
        for(s in 1:length(sample_sizes)){
          norm_control_list[[index]] <- norm.control.lister(trials = nsim, n = sample_sizes[s], 
                                                            p = round((parameters[p] * sample_sizes[s]), 0),
                                                            lxvar = lower_vars[v], uxvar = upper_vars[v], 
                                                            lcor = lower_cors[c], ucor = upper_cors[c],
                                                            sd_in = errors[e])
          index <- index + 1
        }
      }
    }
  }
}

## Running Simulations
norm_output_list <- list()

for(i in 1:simulations){
  norm_output_list[[i]] <- do.call(norm.simulator, norm_control_list[[i]])
}

save(norm_output_list, file = "~/simulations/normal_simulations.rda")

# Non-Normal Simulations ----
## Control vectors
bg_alpha <- c(1, 5)
bg_beta <- c(5, 10)
logistic_vars <- c(1, 10)
t_df <- c(1, 20)

simulations2 <- length(bg_alpha) * length(logistic_vars) * length(t_df) * length(errors) *length(lower_cors) * length(lower_vars)

non_norm_control_list <- list()

index2 <- 1

for(g in 1:length(bg_alpha)){
  for(l in 1:length(logistic_vars)){
    for(df in 1:length(t_df)){
      for(e in 1:length(errors)){
        for(c in 1:length(lower_cors)){
          for(v in 1:length(lower_vars)){
            non_norm_control_list[[index2]] <- 
              non.norm.control.lister(trials = nsim, n = 100, p = 5, lxvar = lower_vars[v], 
                                  uxvar = upper_vars[v], lcor = lower_cors[c], 
                                  ucor = upper_cors[c], sd_in = errors[e],
                                  bgl = bg_alpha[g], bgu = bg_beta[g], 
                                  logis = logistic_vars[l], df = t_df[df])
            index2 <- index2 + 1
          }
        }
      }
    }
  }
}

## Running Simulations
non_norm_output_list <- list()

for(i in 1:simulations2){
  non_norm_output_list[[i]] <- do.call(non.norm.simulator, non_norm_control_list[[i]])
}

save(non_norm_output_list, file = "~/simulations/non_normal_simulations.rda")
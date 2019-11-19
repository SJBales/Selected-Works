library(MASS)
library(propagate)
library(Renvlp)

X.simulator <- function(n, p, lmu, umu, lxvar, uxvar, lcor, ucor) {
  
  # Constructing means and covariance matrix
  mu <- runif(p, lmu, umu)
  cor_vec <- runif(p, sqrt(lcor), sqrt(ucor))
  cor_mat <- cor_vec %*% t(cor_vec)
  x_variances <- runif(p, lxvar, uxvar)
  diag(cor_mat) <- 1
  sigma <- cor2cov(cor_mat, x_variances)

  # Creating X matrix
  X <- mvrnorm(n, mu, sigma, empirical = T)
  return(X)
}

Y.simulator <- function(X, n, sd_in, beta_vec){
  # Creating Responses
  noise <- rnorm(n, mean = 0, sd = sd_in) # increased noise
  Y <- (X %*% unlist(beta_vec)) + noise
  return(Y)
}

sim.ols.estimator <- function(X, Y, p, n, beta_vec) {
  # Linear model
  lin_model <- lm(Y ~ X)
  lin_sum_df <- data.frame(summary(lin_model)$coefficients)

  # Differences from true betas
  lin_beta_diff <- lin_sum_df$Estimate[-1] - beta_vec
  
  # Differences from true betas for OLS
  return(lin_beta_diff)
}

sim.env.estimator <- function(X, Y, p, n, beta_vec) {
  # Envelope Model
  dim_en <- u.xenv(X, Y)
  env_model <- xenv(X, Y, dim_en$u.aic)
  
  # Differences from true betas
  env_beta_diff <- env_model$beta - beta_vec
  
  # Differences from true betas for envelope
  return(c(env_beta_diff))
}

env.ols.simulator <- function(trials, n, p, lbeta, ubeta, lmu, umu, lxvar, uxvar, lcor, ucor, noise){
  
  # Empty Dataframes
  env_diff <- data.frame(matrix(nrow = trials, ncol = p))
  ols_diff <- data.frame(matrix(nrow = trials, ncol = p))
  
  for (i in 1:trials) {
    beta_vec <- runif(p, lbeta, ubeta)
    X <- X.simulator(n, p, lmu, umu, lxvar, uxvar, lcor, ucor)
    Y <- Y.simulator(X, n, noise, beta_vec)
    ols_diff[i, ] <- sim.ols.estimator(X, Y, p, n, beta_vec)
    env_diff[i, ] <- sim.env.estimator(X, Y, p, n, beta_vec)
  }
  
  return(list("OLS" = ols_diff, "Env" = env_diff))
}

sim_control_ols <- list("trials" = 1000,
                        "n" = NA, 
                        "p" = NA, 
                        "lbeta" = 0,
                        "ubeta" = 3,
                        "lmu" = 5, 
                        "umu" = 15, 
                        "lxvar" = 9, 
                        "uxvar" = 225,
                        "lcor" = 0.75,
                        "ucor" = 1,
                        "noise" = 10)

sample_sizes <- c(25, 50, 100, 300)

sim_controls <- rep(list(sim_control_ols), length(sample_sizes))
ols_output <- rep(list(NA), length(sample_sizes))
env_output <- rep(list(NA), length(sample_sizes))
index <- 1

for(j in seq_along(sample_sizes)){
  sim_controls[[index]]$n <- sample_sizes[j]
  sim_controls[[index]]$p <- round((sample_sizes[j] * 0.3), 0)
  ols_output[[index]] <- do.call(env.ols.simulator, sim_controls[[index]])$OLS
  env_output[[index]] <- do.call(env.ols.simulator, sim_controls[[index]])$Env
  index <- index + 1
}

list_output <- list(ols_output, env_output)

save(list_output, file = "~/simulations/ols_sim11.rda")



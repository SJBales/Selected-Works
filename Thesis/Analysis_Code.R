library(tidyverse)
library(Renvlp)
library(ks)
library(xts)
library(ggridges)

load("non_normal_simulations.rda")
load("normal_simulations.rda")

# Simulation Functions----

## Computes MSE
mse.fun <- function(x) {(t(x) %*% x) / length(x)}

## Computed MC SE estimate
mse.mc.se <- function(df_list, x){
  
  unlisted_dfs <- unlist(df_list[x], recursive = F)
  calc <- function(x) {
    mse_est <- mse.fun(x)
    return(sqrt(sum(((x^2) - c(mse_est))^2) / ((1000*999))))
  }
  
  a <- sapply(unlisted_dfs$aic, calc)
  b <- sapply(unlisted_dfs$bic, calc)
  l <- sapply(unlisted_dfs$lrt, calc)
  o <- sapply(unlisted_dfs$ols, calc)
  
  return(vec(cbind(a, b, l, o)))
}

## Converts characters to numeric values
converter <- function(x) {as.numeric(as.character(x))}

## Creates a dataframe summarizing the performance of each method under the simulations scenarios
framer <- function(df_list){
  
  aic_df <- data.frame(cbind(Method = rep("AIC", ncol(df_list$aic)), 
                             Variable = seq(1:ncol(df_list$aic)), 
                             Bias = sapply(df_list$aic, mean),
                             Variance = sapply(df_list$aic, var),
                             MSE = sapply(df_list$aic, mse.fun)),
                       stringsAsFactors = F,
                       row.names = NULL) %>%
    dplyr::mutate_at(.vars = 3:5, .funs = converter)
  
  bic_df <- data.frame(cbind(Method = rep("BIC", ncol(df_list$bic)), 
                             Variable = seq(1:ncol(df_list$bic)), 
                             Bias = sapply(df_list$bic, mean),
                             Variance = sapply(df_list$bic, var),
                             MSE = sapply(df_list$bic, mse.fun)),
                       stringsAsFactors = F,
                       row.names = NULL) %>%
    dplyr::mutate_at(.vars = 3:5, .funs = converter)
  
  lrt_df <- data.frame(cbind(Method = rep("LRT", ncol(df_list$lrt)), 
                             Variable = seq(1:ncol(df_list$lrt)), 
                             Bias = sapply(df_list$lrt, mean),
                             Variance = sapply(df_list$lrt, var),
                             MSE = sapply(df_list$lrt, mse.fun)),
                       stringsAsFactors = F,
                       row.names = NULL) %>%
    dplyr::mutate_at(.vars = 3:5, .funs = converter)
  
  ols_df <- data.frame(cbind(Method = rep("OLS", ncol(df_list$ols)), 
                             Variable = seq(1:ncol(df_list$ols)), 
                             Bias = sapply(df_list$ols, mean),
                             Variance = sapply(df_list$ols, var),
                             MSE = sapply(df_list$ols, mse.fun)),
                       stringsAsFactors = F,
                       row.names = NULL) %>%
    dplyr::mutate_at(.vars = 3:5, .funs = converter)
  
  big_df <- rbind(aic_df, bic_df, lrt_df, ols_df) %>%
    mutate(Method = as.factor(Method))
  
  return(big_df)
}

## Formats dataframes from the output list for import into framer
framing.hammer <- function(dfs, x){
  unlist(dfs[x], recursive = F) %>%
    framer()
}

## Constructs a barplot of MSE for each method (summarized across variables)
mse.summary.plotter <- function(df){
  df %>%
    group_by(Method) %>%
    summarize(MSE_avg = mean(MSE)) %>%
    ggplot(aes(x = Method, y = MSE_avg, fill = Method)) + 
    geom_bar(stat = "identity")
}

# Constructs a barplot of Bias for each method (summarized across variables)
bias.summary.plotter <- function(df){
  df %>%
    group_by(Method) %>%
    summarize(bias_avg = mean(Bias)) %>%
    ggplot(aes(x = Method, y = bias_avg, fill = Method)) + 
    geom_bar(stat = "identity")
}

# Constructs a barplot of Variance for each method (summarized across variables)
variance.summary.plotter <- function(df){
  df %>%
    group_by(Method) %>%
    summarize(variance_avg = mean(Variance)) %>%
    ggplot(aes(x = Method, y = variance_avg, fill = Method)) + 
    geom_bar(stat = "identity")
}

## Computes Monte Carlo Standard Errors
mc.se <- function(df, trial){
  aug_df <- df %>%
    mutate(Bias_MC_SE = sqrt((1 / 1000) * Variance),
           Variance_MC_SE = sqrt(Variance / (2 * (1000 - 1))),
           MSE_MC_SE = mse.mc.se(norm_output_list, trial))
  return(aug_df)
}

## Creating an a dataframe of simulation conditions
ss_vec <- rep(c(rep(25, 4*2), rep(75, 4*8), rep(250, 4*25), rep(25, 4*8), rep(75, 4*22), rep(250, 4*75)), 12)
param_vec <- rep(c(rep(0.1, 4*(2+8+25)), rep(0.3, 4*(8+22+75))), 12)
pred_var_vec <- rep(c(rep("low", 4*((2+8+25) + (8+22+75))), rep("high", 4*((2+8+25)+(8+22+75)))), 6)
correlation_vec <- rep(c(rep("low", 4*2*((2+8+25) + (8+22+75))), rep("medium", 4*2*((2+8+25) + (8+22+75))), rep("high", 4*2*((2+8+25) + (8+22+75)))), 2)

noise_vec <- c(rep('low', (840*4)), rep("high", (840*4)))

simulation_parameters <- data.frame(Sample_size = ss_vec, 
                                    Parameters = param_vec, 
                                    Predictor_Variance = pred_var_vec, 
                                    Correlation = correlation_vec, 
                                    Noise = noise_vec)

## Checking the parameter control dataframe
simulation_parameters %>%
  group_by(Noise, Correlation, Predictor_Variance, Parameters, Sample_size) %>%
  summarize(n())

## Creating a dataframe of aggregate results
aggregate_results <- data.frame()

for(i in 1:72){
  new_df <- framing.hammer(norm_output_list, i)
  aggregate_results <- bind_rows(aggregate_results, new_df)
}

## Binding the simulation control dataframe to the aggregate results dataframe
analyzing_df <- bind_cols(aggregate_results, simulation_parameters) %>%
  mutate(Sample_size = factor(Sample_size),
         Predictor_Variance = factor(Predictor_Variance),
         Parameters = factor(Parameters))

# Fama-Macbeth Functions----

## Step one: time series regressions of asset returns and factors
stage1.fm.estimator <- function(return_df, factor_df, criteria = 1){
  
  beta_env_df <- matrix(nrow = ncol(return_df), ncol = ncol(factor_df))
  ratio_env_df <- matrix(nrow = ncol(return_df), ncol = ncol(factor_df))
  dimension_df <- matrix(nrow = ncol(return_df), ncol = 3) 
  
  X <- data.frame(factor_df, row.names = NULL)
  
  for (s in 1:ncol(return_df)) {
    Y <- data.frame(return_df[, s])
    
    dims <- u.xenv(X, Y)
    
    if (dims[[criteria]] != 0){
      dim <- dims[[criteria]]
    } else {
      dim <- 1
    }
    
    env_model <- xenv(X, Y, u = dim)
    beta_env_df[s,] <- c(env_model$beta)
    ratio_env_df[s,] <- c(env_model$ratio)
    dimension_df[s,] <- c(unlist(dims)[1:3])
  }
  return(list("Beta" = beta_env_df,
              "Ratio" = ratio_env_df,
              "Dimensions" = dimension_df))
}

na.counter <- function(x){sum(is.na(x))}


## Step two: cross-section asset returns as a function of betas estimated in step one
stage2.fm.estimator <- function(beta_df, return_df, criteria = 1){
  
  loadings <- matrix(nrow = nrow(return_df), ncol = ncol(beta_df))
  loadings_ratios <- matrix(nrow = nrow(return_df), ncol = ncol(beta_df))
  dimension_df2 <- matrix(nrow = nrow(return_df), ncol = 3)
  
  for (r in 1:nrow(return_df)) {
    Y2 <- data.frame(t(return_df[r,]), row.names = NULL)
    
    dims2 <- u.xenv(beta_df, Y2)
    
    if (dims2[[criteria]] != 0){
      dim2 <- dims2[[criteria]]
    } else {
      dim2 <- 1
    }
    
    env_model2 <- xenv(beta_df, Y2, u = dim2)
    loadings[r,] <- c(env_model2$beta)
    loadings_ratios[r,] <- c(env_model2$ratio)
    dimension_df2[r,] <- c(unlist(dims2)[1:3])
  }
  return(list("Loadings" = loadings,
              "Loading_Ratios" = loadings_ratios,
              "Dimensions" = dimension_df2))
}

# Normal Simulation Results

# Normal Simulation Plots
analyzing_df %>%
  filter(Sample_size == 250, Predictor_Variance == 'high', Noise == 'low') %>%
  mutate(Correlation = fct_relevel(Correlation, c('low', 'medium', 'high')),
         Parameters = fct_recode(Parameters, '10%' = '0.1', '30%' = '0.3')) %>%
  group_by(Method, Parameters, Correlation) %>%
  summarize(Bias = mean(Bias), SE = mean(Variance), MSE = mean(MSE)) %>%
  ggplot(aes(x = Method, y = MSE, fill = Method)) + 
  facet_wrap(Parameters ~ Correlation) +
  geom_bar(stat = 'identity') +
  ylab("MSE") +
  labs(title = "Figure 1: Large Sample Size",
       subtitle = "High Predictor Variance, Low Response Variation") +
  theme_classic()

analyzing_df %>%
  filter(Sample_size == 75, Predictor_Variance == 'high', Noise == 'low') %>%
  mutate(Correlation = fct_relevel(Correlation, c('low', 'medium', 'high')),
         Parameters = fct_recode(Parameters, '10%' = '0.1', '30%' = '0.3')) %>%
  group_by(Method, Parameters, Correlation) %>%
  summarize(Bias = mean(Bias), SE = mean(Variance), MSE = mean(MSE)) %>%
  ggplot(aes(x = Method, y = MSE, fill = Method)) + 
  facet_wrap(Parameters ~ Correlation) +
  geom_bar(stat = 'identity') +
  ylab("MSE") +
  labs(title = "Figure 2: Medium Sample Size",
       subtitle = "High Predictor Variance, Low Response Variation") +
  theme_classic()

analyzing_df %>%
  filter(Sample_size == 25, Predictor_Variance == 'high', Noise == 'low') %>%
  mutate(Correlation = fct_relevel(Correlation, c('low', 'medium', 'high')),
         Parameters = fct_recode(Parameters, '10%' = '0.1', '30%' = '0.3')) %>%
  group_by(Method, Parameters, Correlation) %>%
  summarize(Bias = mean(Bias), SE = mean(Variance), MSE = mean(MSE)) %>%
  ggplot(aes(x = Method, y = MSE, fill = Method)) + 
  facet_wrap(Parameters ~ Correlation) +
  geom_bar(stat = 'identity') +
  ylab("MSE") +
  labs(title = "Figure 3: Small Sample Size",
       subtitle = "High Predictor Variance, Low Response Variation") +
  theme_classic()

analyzing_df %>%
  filter(Sample_size == 250, Predictor_Variance == 'high', Noise == 'high') %>%
  mutate(Correlation = fct_relevel(Correlation, c('low', 'medium', 'high')),
         Parameters = fct_recode(Parameters, '10%' = '0.1', '30%' = '0.3')) %>%
  group_by(Method, Parameters, Correlation) %>%
  summarize(Bias = mean(Bias), SE = mean(Variance), MSE = mean(MSE)) %>%
  ggplot(aes(x = Method, y = MSE, fill = Method)) + 
  facet_wrap(Parameters ~ Correlation) +
  geom_bar(stat = 'identity') +
  ylab("MSE") +
  labs(title = "Figure 4: Large Sample Size",
       subtitle = "High Predictor Variance, High Response Variation") +
  theme_classic()

analyzing_df %>%
  filter(Sample_size == 25, Predictor_Variance == 'high', Noise == 'high') %>%
  mutate(Correlation = fct_relevel(Correlation, c('low', 'medium', 'high')),
         Parameters = fct_recode(Parameters, '10%' = '0.1', '30%' = '0.3')) %>%
  group_by(Method, Parameters, Correlation) %>%
  summarize(Bias = mean(Bias), SE = mean(Variance), MSE = mean(MSE)) %>%
  ggplot(aes(x = Method, y = MSE, fill = Method)) + 
  facet_wrap(Parameters ~ Correlation) +
  geom_bar(stat = 'identity') +
  ylab("MSE") +
  labs(title = "Figure 5: Small Sample Size",
       subtitle = "High Predictor Variance, High Response Variation") +
  theme_classic()

analyzing_df %>%
  filter(Sample_size == 250, Predictor_Variance == 'low', Noise == 'low') %>%
  mutate(Correlation = fct_relevel(Correlation, c('low', 'medium', 'high')),
         Parameters = fct_recode(Parameters, '10%' = '0.1', '30%' = '0.3')) %>%
  group_by(Method, Parameters, Correlation) %>%
  summarize(Bias = mean(Bias), SE = mean(Variance), MSE = mean(MSE)) %>%
  ggplot(aes(x = Method, y = MSE, fill = Method)) + 
  facet_wrap(Parameters ~ Correlation) +
  geom_bar(stat = 'identity') +
  ylab("MSE") +
  labs(title = "Figure 6: Large Sample Size",
       subtitle = "Low Predictor Variance, Low Response Variation") +
  theme_classic()

analyzing_df %>%
  filter(Sample_size %in% c(25, 75), Predictor_Variance == 'low', Noise == 'low') %>%
  mutate(Correlation = fct_relevel(Correlation, c('low', 'medium', 'high')),
         Sample_size = factor(Sample_size)) %>%
  group_by(Method, Sample_size) %>%
  summarize(Bias = mean(Bias), SE = mean(Variance), MSE = mean(MSE)) %>%
  ggplot(aes(x = Method, y = MSE, fill = Method)) + 
  facet_wrap(. ~ Sample_size, scales = "free_y", nrow = 1) +
  geom_bar(stat = 'identity') +
  ylab("MSE") +
  labs(title = "Figure 7: Small and Medium Sample Size",
       subtitle = "Low Predictor Variance, Low Response Variation") +
  theme_classic()

analyzing_df %>%
  filter(Predictor_Variance == 'low', Noise == 'high') %>%
  mutate(Correlation = fct_relevel(Correlation, c('low', 'medium', 'high')),
         Sample_size = factor(Sample_size)) %>%
  group_by(Method, Sample_size) %>%
  summarize(Bias = mean(Bias), SE = mean(Variance), MSE = mean(MSE)) %>%
  ggplot(aes(x = Method, y = MSE, fill = Method)) + 
  facet_wrap(. ~ Sample_size, scales = "free_y", nrow = 1) +
  geom_bar(stat = 'identity') +
  ylab("MSE") +
  labs(title = "Figure 8: Advantages of Envelopes",
       subtitle = "Low Predictor Variance, High Response Variation") +
  theme_classic()

analyzing_df %>%
  filter(Sample_size == 25, Predictor_Variance == 'low', Noise == 'high') %>%
  mutate(Correlation = fct_relevel(Correlation, c('low', 'medium', 'high')),
         Sample_size = factor(Sample_size)) %>%
  group_by(Method, Sample_size) %>%
  summarize(Bias = mean(Bias), SE = mean(Variance), MSE = mean(MSE))

analyzing_df %>%
  filter(Sample_size == 75, Parameters == 0.1, Predictor_Variance == 'low', 
         Correlation == 'high', Noise == 'low') %>%
  group_by(Method) %>%
  summarize(Bias = mean(Bias), SE = mean(Variance), MSE = mean(MSE))  %>%
  ggplot(aes(x = Method, y = MSE, fill = Method)) +
  geom_bar(stat = 'identity') +
  labs(title = "Figure 9: Advantages of BIC", 
       subtitle = "75 Observations, 8 Highly Collinear Predictors with Low Variability, Low Variability in the Response") +
  theme_classic()

# Non-Normal Simulation Results----
# Creating a dataframe of simulation results
non_normal_aggregate_df <- data.frame()

for(j in 1:96){
  nn_new_df <- framing.hammer(non_norm_output_list, j)
  non_normal_aggregate_df <- bind_rows(non_normal_aggregate_df, nn_new_df)
}

# Creating a dataframe of simuation parameters
norm_var <- rep(c(rep("low", 9*4), rep('high', 9*4)), 48)
norm_cor <- rep(c(rep("low", 2*9*4), rep("medium", 2*9*4), rep("high", 2*9*4)), 16)
norm_error <- rep(c(rep("low", 2*9*4*2), rep("high", 2*9*4*2)), 12)
t_df <- rep(c(rep("1", 2*9*4*2*2), rep("20", 2*9*4*2*2)), 6)
logistic_var <- rep(c(rep("low", 2*9*4*2*2*2), rep("high", 2*9*4*2*2*2)), 3)
beta_var <- rep(c(rep("low", 2*9*4*2*2*2*3), rep("high", 2*9*4*2*2*2*3)))

non_normal_simulation_parameters <- data.frame(norm_var, norm_cor, norm_error, t_df, logistic_var, beta_var)

# Creating an analyzing dataframe
nn_analyzing_df <- data.frame(non_normal_aggregate_df, non_normal_simulation_parameters)

# Overall MSE of each selection criteria
nn_analyzing_df %>% 
  group_by(Method) %>%
  summarize(MSE = sum(MSE), Bias = sum(Bias), Variance = sum(Variance)) %>%
  ggplot(aes(x = Method, y = MSE, fill = Method)) + 
  geom_bar(stat = "identity")

# MSE of Selection Criteria Across Variables
nn_analyzing_df %>%
  group_by(Method, Variable) %>%
  summarize(MSE = mean(MSE)) %>%
  mutate(Variable = converter(Variable)) %>% 
  ggplot(aes(x = Variable, y = MSE, color = Method)) +
  geom_line()

# MSE of Selection Criteria Across vars except gamma and beta
nn_analyzing_df %>%
  group_by(Method, Variable) %>%
  filter(Variable != 6 & Variable != 7) %>%
  summarize(MSE = mean(MSE)) %>%
  ggplot(aes(x = Variable, y = MSE, fill = Method)) +
  geom_bar(stat = 'identity', position = 'dodge')

# MSE of criteria when the variability of normal predictors is high
nn_analyzing_df %>%
  group_by(Method, Variable) %>%
  filter(Variable != 6 & Variable != 7, 
         norm_var == 'high') %>%
  summarize(MSE = mean(MSE)) %>%
  ggplot(aes(x = Variable, y = MSE, fill = Method)) +
  geom_bar(stat = 'identity', position = 'dodge')

# MSE Ridgeline for 
nn_analyzing_df %>%
  filter(Variable == 7) %>%
  ggplot(aes(x = MSE, y = Method, fill = Method)) +
  geom_density_ridges()

# Histograms of MSE for T dist
nn_analyzing_df %>%
  filter(Variable == 7) %>%
  ggplot(aes(x = MSE, fill = Method)) +
  geom_histogram() +
  facet_wrap(Method ~. )

# Application----

# Reading and formatting data
returns <- read.csv("sp_returns.csv") %>%
  mutate(Date = as.Date(as.character(Date), format = "%Y - %m - %d")) %>%
  select_if(colSums(is.na(.))==0)

factors <- read.csv("five_factors.csv") %>%
  mutate(Date = as.Date(as.character(Date), format = "%Y - %m - %d"))

ret_xts <- as.xts(returns[,-1], order.by = returns[,1])
fact_xts <- as.xts(factors[,-1], order.by = factors[,1])

periodicity(ret_xts["/2019-11-01"])
periodicity(fact_xts["1987-06-30/"])

ret_reg <- ret_xts["/2019-11-01"]
fact_reg <- fact_xts["1987-06-30/"]

# Predictor Variability
fact_reg %>%
  data.frame(row.names = NULL) %>%
  summarize_all(.funs = sd) %>%
  t()

# ANOVA
aapl <- data.frame(ret_reg[, 1], fact_reg, row.names = NULL)
aapl_lm <- lm(AAPL ~., data = aapl)
sd(aapl_lm$residuals)

# Correlation Matrix
fact_reg %>% cor()

# Standard Deviation of Betas--Used as predictors in second stage
beta_df %>%
  summarize_all(.funs = sd) %>%
  t()

# Correlation of Betas
beta_df %>%
  cor()

# MSE from stage 2 calc with betas
stage2_summary_df <- data.frame(t(ret_reg[1,]), beta_df, row.names = NULL) %>%
  rename(Return = 1)

sd(lm(Return ~., data = stage2_summary_df)$residuals)

sd(stage2_summary_df$Return)

# AIC
stage1_results <- stage1.fm.estimator(ret_reg, fact_reg, criteria = 1)
beta_df <- data.frame(stage1_results$Beta)
stage2_results <- stage2.fm.estimator(beta_df, ret_reg, criteria = 1)

stage2_results$Loading_Ratios %>%
  data.frame() %>%
  summarize_all(.funs = mean)

# BIC
stage1_results_bic <- stage1.fm.estimator(ret_reg, fact_reg, criteria = 2)
beta_df_bic <- data.frame(stage1_results_bic$Beta)
stage2_results_bic <- stage2.fm.estimator(beta_df_bic, ret_reg, criteria = 2)

stage2_results_bic$Loading_Ratios %>%
  data.frame() %>%
  summarize_all(.funs = mean)

# LRT
stage1_results_lrt <- stage1.fm.estimator(ret_reg, fact_reg, criteria = 3)
beta_df_lrt <- data.frame(stage1_results_lrt$Beta)
stage2_results_lrt <- stage2.fm.estimator(beta_df_bic, ret_reg, criteria = 3)

stage2_results_lrt$Loading_Ratios %>%
  data.frame() %>%
  summarize_all(.funs = mean)
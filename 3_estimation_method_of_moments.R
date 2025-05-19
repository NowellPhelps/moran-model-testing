rm(list = ls())
library(tidyverse)

load("biased_1000_sims_w_55.RData")

n_individuals <- 100

get_alpha_hat_k <- function(mean_change, prop, n_individuals){
  f_a <- prop
  f_b <- 1- prop
  
  num <- mean_change/f_b + 1
  denom <- num + 1 - mean_change/f_a
  
  alpha_hat_k <- num/denom
  # 
  # alpha_hat_k[which(alpha_hat_k > 1)] <- 1
  # alpha_hat_k[which(alpha_hat_k <0 )] <- 0
  # 
  
  return(alpha_hat_k)
}

get_alpha_hat_estimator <- function(trajectory, n_individuals, my_sim_id, bootstrap = F){ 
  
  sim <- trajectory 
  if(!(bootstrap)){
    sim <- sim %>% 
      filter(sim_id == my_sim_id) %>%
      mutate(lead = lead(prop), 
             change = 100*(lead - prop))
    
    max_it <- max(sim$it)
    
    sim <- sim %>%
      filter(it < max_it) 
  }
  
  sim <- sim %>%
    group_by(prop) %>%
    summarise(mean_change = mean(change), 
              n = n()) %>%
    filter(!(prop %in% c(0, 1)))
  
  n_total <- sum(sim$n)
  sim <- sim
  sim$alpha_hat <- get_alpha_hat_k(sim$mean_change, sim$prop, n_individuals)
  
  alpha_hat <- sum(sim$n*sim$alpha_hat)/n_total
  return(alpha_hat)
}

get_one_bootstrap <- function(sim, n_individuals){
  sim_sample <- sim[sample(1:nrow(sim), nrow(sim), replace = T),]
  return(get_alpha_hat_estimator(sim_sample, n_individuals, bootstrap = T))
}


bootstrap_alpha_hat_estimator <- function(trajectory, n_bootstraps, n_individuals, my_sim_id){
  
  sim <- trajectory %>% 
    filter(sim_id == my_sim_id) %>%
    mutate(lead = lead(prop), 
           change = 100*(lead - prop))
  max_it <- max(sim$it)
  sim <- sim %>%
    filter(it < max_it)
  
  bootstrapped_estimates <- replicate(n_bootstraps, get_one_bootstrap(sim, 100))
  
  return(bootstrapped_estimates)
}

alphas <- sapply(1:1000, FUN = function(x) get_alpha_hat_estimator(sims, 100, x))

bootstrap_coverage <- function(sims, n_bootstraps = 100, n_individuals = 100, my_sim_id, true_alpha = 0.55){
  bootstrap_test <- bootstrap_alpha_hat_estimator(sims, n_bootstraps = 1000, 100, my_sim_id)
  coverage <- ifelse(quantile(bootstrap_test, 0.025) < true_alpha & 
                       quantile(bootstrap_test, 0.975) > true_alpha, 1, 0)
  return(coverage)
}

bootstrap_test <- function(sims, n_bootstraps = 100, n_individuals = 100, my_sim_id, true_alpha = 0.55){
  bootstrap_test <- bootstrap_alpha_hat_estimator(sims, n_bootstraps = 1000, 100, my_sim_id)
  reject <- ifelse(quantile(bootstrap_test, 0.025) > 0.05, 1, 0)
  return(reject)
}

ptm <- proc.time()
coverage <- parallel::mclapply(1:1000, FUN = function(x) bootstrap_coverage(sims, 100, 100, x, true_alpha = 0.55))
rejections <- sapply(1:10, FUN = function(x) bootstrap_test(sims, 100, 100, x, true_alpha = 0.55))
print(ptm - proc.time())







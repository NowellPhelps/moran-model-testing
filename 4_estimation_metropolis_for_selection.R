rm(list = ls())
library(tidyverse)
library(grid)
library(gridExtra)
source("allele_functions.R")

load("unbiased_1000_sims.RData")

n_individuals <- 100

get_preprocessed_trajectory <- function(sims, my_sim_id){
  sim <- sims %>% 
    filter(sim_id == my_sim_id) %>%
    mutate(lead = lead(prop), 
           change = round(100*(lead - prop)))
  max_it <- max(sim$it)
  sim <- sim %>%
    filter(it < max_it) %>%
    group_by(prop) %>%
    summarise(n_up = length(which(change  == 1)),
              n_down = length(which(change == -1)),
              n_stay = length(which(change == 0)),
              n = n()) %>%
    filter(!(prop %in% c(0, 1)))
  
  return(sim)
}

get_lhood_one_row_trajectory <- function(prop, n_up, n_down, n_stay, transition_matrix){
  state_id <- round(prop*100) + 1
  return(transition_matrix[state_id, state_id]^n_stay*transition_matrix[state_id, state_id-1]^n_down*transition_matrix[state_id, state_id+1]^n_up)
}

calculate_log_likelihood <- function(alpha_current, sim, n_individual){
  transition_matrix <- get_transition_matrix(n_individual, alpha_current)
  lhood <- 0
  for(i in 1:nrow(sim)){
    prop <- sim$prop[i]
    n_up <- sim$n_up[i]
    n_down <- sim$n_down[i]
    n_stay <- sim$n_stay[i]
    lhood <- lhood + log(get_lhood_one_row_trajectory(prop, n_up, n_down, n_stay, transition_matrix))
  }
  return(lhood)
}



run_one_chain <- function(maxIts = 1000, initial_alpha = 0.4, my_sim_id, sims){
  
  sim <- get_preprocessed_trajectory(sims, my_sim_id)
  
  mcmc_traj <- c()
  alpha_current <- initial_alpha
  log_lhood_current <- calculate_log_likelihood(alpha_current, sim, 100)
  mcmc_traj <- c(alpha_current)
  
  
  for(i in 1:maxIts){
    alpha_prop <- rnorm(1, alpha_current, 0.08)
    if (alpha_prop > 1 | alpha_prop < 0){
      alpha_prop <- alpha_current
    } else{
      log_lhood_proposed <- calculate_log_likelihood(alpha_prop, sim, 100)
      if(log_lhood_proposed > log_lhood_current){
        alpha_current <- alpha_prop
        log_lhood_current <- log_lhood_proposed
      } else{
        if (runif(1) < exp(log_lhood_proposed - log_lhood_current)){
          alpha_current <- alpha_prop
          log_lhood_current <- log_lhood_proposed
        } 
      }
      mcmc_traj <- c(mcmc_traj, alpha_current)
    }
  }
  return(mcmc_traj)
}


make_traceplots <- function(my_sim_id, trajectories, true_alpha){
  
  plot_data <- data.frame(alpha = trajectories[[(my_sim_id)]], 
                          it = c(0:(length(trajectories[[(my_sim_id)]])-1)))
  
  traceplot <- ggplot(data = plot_data, aes(x = it, y = alpha, colour = ifelse(it >= 100, "Keep", "Discard"))) + 
    geom_line() +
    scale_colour_discrete(name = "Iteration") + 
    labs(y = "Alpha", 
         x = "Iteration") +
    scale_y_continuous(limits = c(0.3, 0.8)) +
    geom_hline(yintercept = true_alpha, colour = "red") +
    geom_hline(yintercept = mean(plot_data$alpha[which(plot_data$it >= 100)]), linetype = "dashed", colour = "grey20") +
    geom_hline(yintercept = quantile(plot_data$alpha[which(plot_data$it >= 100)], 0.025), linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = quantile(plot_data$alpha[which(plot_data$it >= 100)], 0.975), linetype = "dashed", colour = "grey50") +
    theme_classic()
  
  return(traceplot)
}

get_coverage <- function(my_sim_id, trajectories, true_alpha){
  plot_data <- data.frame(alpha = trajectories[[(my_sim_id)]], 
                          it = c(0:(length(trajectories[[(my_sim_id)]])-1)))
  
  q_025 <- quantile(plot_data$alpha[which(plot_data$it >= 100)], 0.025)
  q_975 <- quantile(plot_data$alpha[which(plot_data$it >= 100)], 0.975)
  
  coverage <- ifelse(true_alpha <= q_975 & true_alpha >= q_025, 1, 0)
  return(coverage)
}

trajectories <- lapply(1:1000, FUN = function(my_sim_id) run_one_chain(maxIts = 2000, initial_alpha = 0.4, my_sim_id, sims))
save(trajectories, file = "mcmc_chains_50.RData")
ps <- lapply(1:1000, FUN = function(my_sim_id) make_traceplots(my_sim_id, trajectories, true_alpha = 0.5))


save(trajectories, file = "mcmc_chains_50.RData")

load("mcmc_chains_55.RData")
cairo_pdf(paste0("traceplots_05.pdf"), height = 5, width = 10, onefile=T)

for(my_sim_id in 1:1000){
  grid.arrange(ps[[my_sim_id]])
}
dev.off()


coverage_2 <- lapply(1:1000, FUN = function(my_sim_id) get_coverage(my_sim_id, trajectories, true_alpha = 0.55))
table(unlist(coverage_2))



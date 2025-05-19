rm(list = ls())
library(tidyverse)
library(grid)
library(gridExtra)

source("many_genes_functions.R")

traj <- run_one_trajectory_many_genes(nIts = 10000, n_individuals = 100, n_genes = 2, weights = c(0.8, 0.2), 
                                      reproduction_type = "uniform", rho = .999, seed = NULL)


proportions <- traj$trajectory
colnames(proportions) <- gsub("prop.", "", colnames(proportions))
proportions <- proportions %>%
  pivot_longer(cols = as.character(1:2), names_to = c("gene"), values_to = "prop")

p <- ggplot(proportions %>% filter(it < 10000), aes(x = it, y = prop*100, group = gene, colour = gene)) +
  geom_line() + 
  theme_classic() +
  labs(x = "Iteration",
       y = "Number of individuals with Allele A") +
  scale_colour_discrete(name = "Gene")

proportions <- traj$trajectory_proportions_complete

name_ref <- traj$possible_phenotypes %>%
  mutate(name_long <- paste(Var1, Var2))

colnames(proportions) <- gsub("proportions_complete.", "", colnames(proportions))
colnames(proportions) <- c("it", name_ref$name_long)

proportions <- proportions %>%
  pivot_longer(cols = name_ref$name_long, names_to = c("gene"), values_to = "prop")

p <- ggplot(proportions %>% filter(it < 3000), aes(x = it, y = prop*100, group = gene, colour = gene)) +
  geom_line() + 
  theme_classic() +
  labs(y = "Number of individuals", 
       x = "Iteration") +
  scale_colour_discrete(name = "Phenotype")




sim <- traj$trajectory_proportions_complete
names(sim) <- c("it", "prop_BB", "prop_AB", "prop_BA", "prop_AA")

sim  <- sim %>%
  mutate(lead_AA = lead(prop_AA),
         change_AA  = round(100*(lead_AA - prop_AA)),
         lead_AB = lead(prop_AB), 
         change_AB  = round(100*(lead_AB - prop_AB)),
         lead_BA = lead(prop_BA), 
         change_BA  = round(100*(lead_BA - prop_BA)),
         lead_BB = lead(prop_BB), 
         change_BB  = round(100*(lead_BB - prop_BB)))

maxIt <- max(sim$it)

sim <- sim %>% 
  filter(it < maxIt) 

sim <- sim[which(pmax(sim$prop_AA, sim$prop_AB, sim$prop_BA,sim$prop_BB) < 1),]

get_lhood_one_row_trajectory <- function(prop_AA, prop_AB, prop_BA, prop_BB, 
                                         change_AA, change_AB, change_BA, change_BB, rho_current, alpha_1_current, alpha_2_current, n_individuals){
  # Rewrite quantities
  X_AA <- prop_AA*n_individuals
  X_AB <- prop_AB*n_individuals
  X_BA <- prop_BA*n_individuals
  X_BB <- prop_BB*n_individuals

  X_plusB <- X_AB + X_BB
  X_plusA <- X_AA + X_BA
  
  # Probability of reproduction
  p_aa_reproduce <- 1/(n_individuals*(n_individuals-1))*(X_AA*(X_plusA - 1 + rho_current*X_plusB) + (1-rho_current)*X_AB*X_plusA)
  p_ab_reproduce <- 1/(n_individuals*(n_individuals-1))*(X_AB*(X_plusB - 1 + rho_current*X_plusA) + (1-rho_current)*X_AA*X_plusB)
  p_ba_reproduce <- 1/(n_individuals*(n_individuals-1))*(X_BA*(X_plusA - 1 + rho_current*X_plusB) + (1-rho_current)*X_BB*X_plusA)
  p_bb_reproduce <- 1/(n_individuals*(n_individuals-1))*(X_BB*(X_plusB - 1 + rho_current*X_plusA) + (1-rho_current)*X_BA*X_plusB)
  
  # Probability of death
  p_aa_die <- X_AA*(1-alpha_1_current)*(1-alpha_2_current)
  p_ab_die <- X_AB*(1-alpha_1_current)*(alpha_2_current)
  p_ba_die <- X_BA*(alpha_1_current)*(1-alpha_2_current)
  p_bb_die <- X_BB*(alpha_1_current)*(alpha_2_current)
  S <- p_aa_die + p_ab_die + p_ba_die + p_bb_die
  p_aa_die <- p_aa_die/S
  p_ab_die <- p_ab_die/S
  p_ba_die <- p_ba_die/S
  p_bb_die <- p_bb_die/S

  # Calculate row likelihood using data
  if(max(abs(c(change_AA, change_AB, change_BA, change_BB))) == 0){
    # No change - same type of individual born and dies
    p <- p_aa_reproduce*p_aa_die + p_ab_reproduce*p_ab_die + p_ba_reproduce*p_ba_die + p_bb_reproduce*p_bb_die 
  } else{
    p <- 1
    # Change - known phenotype dies and known phenotype is born
    for (i in 1:4){
      if(c(change_AA, change_AB, change_BA, change_BB)[i] == -1){
        p <- p *c(p_aa_die, p_ab_die, p_ba_die, p_bb_die)[i]
      } else if(c(change_AA, change_AB, change_BA, change_BB)[i] == 1){
        p <- p *c(p_aa_reproduce, p_ab_reproduce, p_ba_reproduce, p_bb_reproduce)[i]
      }
    }
  }
  
  return(p)
}


calculate_log_likelihood <- function(rho_current, alpha_1_current, alpha_2_current, sim, n_individual){
  lhood <- 0
  for(i in 1:nrow(sim)){
    lhood <- lhood + log(get_lhood_one_row_trajectory(prop_AA = sim$prop_AA[i],
                                                      prop_AB = sim$prop_AB[i],
                                                      prop_BA = sim$prop_BA[i],
                                                      prop_BB = sim$prop_BB[i],
                                                      change_AA = sim$change_AA[i],
                                                      change_AB = sim$change_AB[i],
                                                      change_BA = sim$change_BA[i],
                                                      change_BB = sim$change_BB[i],
                                                      rho_current = rho_current,
                                                      alpha_1_current = alpha_1_current,
                                                      alpha_2_current = alpha_2_current,
                                                      n_individuals = n_individual))
  }
  return(lhood)
}



run_one_chain <- function(maxIts = 1000, initial_rho = 0.5, initial_alpha_1 = 0.5, initial_alpha_2 = 0.5, sim, prop_sd_rho = NULL, prop_sd_alpha = NULL){
  
  mcmc_traj <- c()
  rho_current <- initial_rho
  alpha_1_current <- initial_alpha_1
  alpha_2_current <- initial_alpha_2
  log_lhood_current <- calculate_log_likelihood(rho_current, alpha_1_current, alpha_2_current, sim, 100)
  mcmc_traj <- data.frame(rho = rho_current,
                          alpha_1 = alpha_1_current,
                          alpha_2 = alpha_2_current)
  
  if(is.null(prop_sd_rho)){
    prop_sd_rho <- 0.1
  } 
  
  if(is.null(prop_sd_alpha)){
    prop_sd_alpha <- 0.04
  }
  
  for(i in 1:maxIts){
    print(i)
    # Update rho
    rho_prop <- 2
    while (rho_prop  > 1 | rho_prop  < 0){
      rho_prop <- rnorm(1, rho_current, prop_sd_rho)
    } 
    
    log_lhood_proposed <- calculate_log_likelihood(rho_prop, alpha_1_current, alpha_2_current, sim, 100)
    if(log_lhood_proposed > log_lhood_current){
      rho_current <- rho_prop
      log_lhood_current <- log_lhood_proposed
    } else{
      if (runif(1) < exp(log_lhood_proposed - log_lhood_current)){
        rho_current <- rho_prop
        log_lhood_current <- log_lhood_proposed
      }
    }
    
    # Update alpha 1 
    alpha_1_prop <- 2
    while (alpha_1_prop  > 1 | alpha_1_prop  < 0){
      alpha_1_prop <- rnorm(1, alpha_1_current, prop_sd_alpha)
    } 
    
    log_lhood_proposed <- calculate_log_likelihood(rho_current, alpha_1_prop, alpha_2_current, sim, 100)
    if(log_lhood_proposed > log_lhood_current){
      alpha_1_current <- alpha_1_prop
      log_lhood_current <- log_lhood_proposed
    } else{
      if (runif(1) < exp(log_lhood_proposed - log_lhood_current)){
        alpha_1_current <- alpha_1_prop
        log_lhood_current <- log_lhood_proposed
      }
    }
    
    # Update alpha 2
    alpha_2_prop <- 2
    while (alpha_2_prop  > 1 | alpha_2_prop  < 0){
      alpha_2_prop <- rnorm(1, alpha_2_current, prop_sd_alpha)
    } 
    
    log_lhood_proposed <- calculate_log_likelihood(rho_current, alpha_1_current, alpha_2_prop, sim, 100)
    if(log_lhood_proposed > log_lhood_current){
      alpha_2_current <- alpha_2_prop
      log_lhood_current <- log_lhood_proposed
    } else{
      if (runif(1) < exp(log_lhood_proposed - log_lhood_current)){
        alpha_2_current <- alpha_2_prop
        log_lhood_current <- log_lhood_proposed
      }
    }
    
    mcmc_traj <- rbind(mcmc_traj, 
                       data.frame(rho = rho_current,
                                  alpha_1 = alpha_1_current,
                                  alpha_2 = alpha_2_current))
    
  }
  return(mcmc_traj)
}

chain <- run_one_chain(maxIts = 5000, initial_rho = 0.5, initial_alpha_1 = 0.5, initial_alpha_2 = 0.5, sim, prop_sd_alpha = 0.1)

true_alpha_1 <- 0.8
true_alpha_2 <- 0.2
true_rho = 0.999
burn_in <- 200

plot_data <- cbind(chain, 0:(nrow(chain)-1))
names(plot_data)[4] <- "it"

traceplot_rho <- ggplot(data = plot_data, aes(x = it, y = rho, colour = ifelse(it >= burn_in, "Keep", "Discard"))) + 
  geom_line() +
  scale_colour_discrete(name = "Iteration") + 
  labs(y = "Rho", 
       x = "Iteration") +
  scale_y_continuous(limits = c(.97, 1.005)) +
  geom_hline(yintercept = true_rho, colour = "purple") +
  geom_hline(yintercept = mean(plot_data$rho[which(plot_data$it >= burn_in)]), linetype = "dashed", colour = "grey20") +
  geom_hline(yintercept = quantile(plot_data$rho[which(plot_data$it >= burn_in)], 0.025), linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = quantile(plot_data$rho[which(plot_data$it >= burn_in)], 0.975), linetype = "dashed", colour = "grey50") +
  theme_classic()

plot(traceplot_rho)


traceplot_alpha_1<- ggplot(data = plot_data, aes(x = it, y = alpha_1, colour = ifelse(it >= burn_in, "Keep", "Discard"))) + 
  geom_line() +
  scale_colour_discrete(name = "Iteration") + 
  labs(y = "Alpha 1", 
       x = "Iteration") +
  scale_y_continuous(limits = c(0, 1)) +
  geom_hline(yintercept = true_alpha_1, colour = "purple") +
  geom_hline(yintercept = mean(plot_data$alpha_1[which(plot_data$it >= burn_in)]), linetype = "dashed", colour = "grey20") +
  geom_hline(yintercept = quantile(plot_data$alpha_1[which(plot_data$it >= burn_in)], 0.025), linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = quantile(plot_data$alpha_1[which(plot_data$it >= burn_in)], 0.975), linetype = "dashed", colour = "grey50") +
  theme_classic()

plot(traceplot_alpha_1)


traceplot_alpha_2<- ggplot(data = plot_data, aes(x = it, y = alpha_2, colour = ifelse(it >= burn_in, "Keep", "Discard"))) + 
  geom_line() +
  scale_colour_discrete(name = "Iteration") + 
  labs(y = "Alpha 2", 
       x = "Iteration") +
  scale_y_continuous(limits = c(0, 1)) +
  geom_hline(yintercept = true_alpha_2, colour = "purple") +
  geom_hline(yintercept = mean(plot_data$alpha_2[which(plot_data$it >= burn_in)]), linetype = "dashed", colour = "grey20") +
  geom_hline(yintercept = quantile(plot_data$alpha_2[which(plot_data$it >= burn_in)], 0.025), linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = quantile(plot_data$alpha_2[which(plot_data$it >= burn_in)], 0.975), linetype = "dashed", colour = "grey50") +
  theme_classic()

plot(traceplot_alpha_2)







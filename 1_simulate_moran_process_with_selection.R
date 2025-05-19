library(tidyverse)
library(doParallel)
library(foreach)
library(dplyr)
library(wesanderson)

allele_fill_scale <- scale_fill_manual(values = c("A" = wes_palette("Darjeeling1")[1],
                                                  "B" = wes_palette("Darjeeling1")[2]))


run_one_trajectory_one_gene <- function(nIts, n_individuals = 100, seed = NULL, initial_pop = NULL, weight = 0.5){
  if (!(is.null(seed))){
    set.seed(seed)
  }
  
  # Initialize population
  if(is.null(initial_pop)){
    pop <- sample(rep(c(0,1), n_individuals/2), n_individuals)
  } else{
    pop <- initial_pop
  }
  
  proportions <- length(which(pop == 1))/length(pop)
  
  # Evolution
  trajectory <- data.frame(it = 0, prop = t(proportions))
  
  for(i in 1:nIts){
    # select individual to die and individual to reproduce
    sample_weights <- ifelse(pop == 1, 1-weight, weight)
    sample_die <- sample(1:n_individuals, 1, prob = sample_weights)
    sample_reproduce <- sample(1:n_individuals, 1) 
    pop <- c(pop, pop[sample_reproduce])
    pop <- pop[-sample_die]
    proportions <- length(which(pop == 1))/length(pop)
    trajectory <- rbind(trajectory, data.frame(it = i, prop = t(proportions)))
  }
  return(trajectory)
}

run_simulations <- function(nSims, nIts = 1000, n_individuals = 100, initial_pop = NULL, weight = 0.5){
  
  simulations <- foreach(sim = 1:nSims, .combine = bind_rows, .packages = c("dplyr")) %dopar% {
    source("allele_functions.R")
    trajectory <- run_one_trajectory_one_gene(nIts, n_individuals, seed = sim, initial_pop, weight) %>%
      mutate(sim_id = sim)
    return(trajectory)
  }
  return(simulations)
}

get_transition_matrix <- function(n_individuals, weight){
  transition_matrix <- NULL
  for (i in c(0, 1:(n_individuals))){
    if(i == 0){
      row <- c(1, rep(0, n_individuals))
    } else if (i == n_individuals){
      row <- c(rep(0, n_individuals), 1)
    } else{
      row <- rep(0, n_individuals +1)
      denom <- (n_individuals - 1)*(weight*(n_individuals - i) + i*(1 - weight))
      row[i+2] <- i*(n_individuals - i)*weight/denom
      row[i+1] <- (weight*(n_individuals - i)*(n_individuals - i -1) + (1-weight)*i*(i-1))/denom
      row[i] <- i*(n_individuals - i)*(1-weight)/denom
    }
    transition_matrix <- rbind(transition_matrix, row)
  }
  return(transition_matrix)
}

get_state_at_time_n <- function(transition_matrix, initial_state, n){
  p <- eigen(transition_matrix)$vectors
  d <- diag(eigen(transition_matrix)$values)
  if (max(abs(transition_matrix - p %*% d %*% solve(p))) > 10e-10){
    stop("Numerical eigendecomposition appears poor")
  }
  return(t(initial_state) %*% p %*% diag((diag(d)**n))%*% solve(p))
}





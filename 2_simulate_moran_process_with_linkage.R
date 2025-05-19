library(MASS)
sample_individual_indicators <- function(n_genes = 2, rho){
  ids <- rep(0, n_genes)
  if(n_genes == 2){
    ids[1] <- rbinom(1,1, 0.5) + 1
    ids[2] <- ifelse(rbinom(1,1, rho) == 1, ids[1], 3 - ids[1])
  } else{
    stop("This number of genes is not yet supported")
  }
  return(ids)
}

reproduce <- function(pop, id1, id2, n_genes = 2, rho){
  ids <- sample_individual_indicators(n_genes, rho)
  child <- rep(0, n_genes)
  child[which(ids == 1)] <- pop[id1, which(ids == 1)]
  child[which(ids == 2)] <- pop[id2, which(ids == 2)]
  
  return(child)
}

run_one_trajectory_many_genes <- function(nIts, n_individuals = 100, n_genes = 2, weights = NULL, 
                                            reproduction_type = "uniform", rho = 0.5, seed = NULL){
  if (!(is.null(seed))){
    set.seed(seed)
  }
  possible_phenotypes <- expand.grid(rep(list(c(0, 1)), n_genes))
  
  # Initialize population
  pop1 <- sample(rep(c(0,1), n_individuals/2), n_individuals)
  pop <- cbind(pop1, pop1)
  colnames(pop) <- c("1", "2")
  
  proportions <- apply(pop, 2, FUN = function(x) length(which(x == 1))/length(x))
  
  case_indices <- apply(pop, 1, FUN = function(x) which(apply(possible_phenotypes, 1, function(row) all(row == x))))
  proportions_complete <- sapply(1:nrow(possible_phenotypes), FUN = function(x) length(which(case_indices == x))/n_individuals)
  
  # Evolution
  trajectory <- data.frame(it = 0, prop = t(proportions))
  trajectory_proportions_complete <-  data.frame(it = 0, proportions_complete = t(proportions_complete))
  
  for(i in 1:nIts){
    
    if (reproduction_type == "uniform"){
      parent_ids <- sample(1:n_individuals, 2, replace = F)
      parent_ids <-  sample(parent_ids, 2, replace = F)
    }
    
    child <- reproduce(pop, parent_ids[1], parent_ids[2], n_genes = 2, rho)
    
    # select individual to die
    if (is.null(weights)){
      sample <- sample(1:n_individuals, 1)
    } else{
      probs <- (pop[,1]*(1-2*weights[1]) + weights[1]) * (pop[,2]*(1-2*weights[2]) + weights[2])
      sample <- sample(1:n_individuals, 1, prob = probs)
    }
    pop <- pop[-sample,]
    
    # select individuals to reproduce
    pop <- rbind(pop, child)
    
    proportions <- apply(pop, 2, FUN = function(x) length(which(x == 1))/length(x))
    
    case_indices <- apply(pop, 1, FUN = function(x) which(apply(possible_phenotypes, 1, function(row) all(row == x))))
    proportions_complete <- sapply(1:nrow(possible_phenotypes), FUN = function(x) length(which(case_indices == x))/n_individuals)
    
    
    trajectory <- rbind(trajectory, data.frame(it = i, prop = t(proportions)))
    trajectory_proportions_complete <- rbind(trajectory_proportions_complete, data.frame(it = i, proportions_complete = t(proportions_complete)))
  }
  
  return(list(trajectory = trajectory,
              possible_phenotypes  =  possible_phenotypes , 
              trajectory_proportions_complete = trajectory_proportions_complete)) 
  
}


run_simulations_multiple_genes <- function(nSims, 
                            nIts = 10000, 
                            n_individuals = 100, weights = NULL, reproduction_type = "uniform",rho = 0.5){
  
  simulations <- foreach(sim = 1:nSims, .combine = bind_rows, .packages = c("dplyr")) %dopar% {
    source("many_genes_functions.R")
    trajectory <- run_one_trajectory_many_genes(nIts, seed = sim, weights = weights, reproduction_type = reproduction_type, rho = rho)$trajectory %>%
      mutate(sim_id = sim)
    return(trajectory)
  }
  return(simulations)
}

  


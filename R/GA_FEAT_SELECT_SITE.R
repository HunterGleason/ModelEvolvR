#'Initializes a population of random subset of sites (rows).
#'
#' This function initializes a population of random sites
#' provided a data frame of possible site locations with 
#' covariates as columns. 
#'
#' @param samp_df Data frame containing all covarites of interest as columns, and all potential sites as rows. 
#' @param num_sites Maximum number of sites to include in a subset, i.e., total number of sample sites.
#' @param N Total number of random subsets of sites to initialize,i.e., population size.
#' @return A list 'population' of random subsets of sites (rows) from the samp_df data.frame. 
#' @export
init_pop_site<-function(N,num_sites,samp_df)
{
  pop<-list()
  for(e in c(1:N))
  {
    pop[[e]]<-sample(c(1:nrow(samp_df)),num_sites)
  }
  
  return(pop)
}

#'Function for calculating 'fitness' of subset of sites based on cumulative kolmogorov smirnov D statistic
#'
#' @param sample Integer vector representing subset of sites (rows) in the candidate site data.frame.
#' @param pop_df Data.frame with same covariates (columns) as the candidate site data.frame, but assumed to represent the population distribution. 
#' @return Cumulative kolmogorov smirnov D statistic for each covaraite in the sample of sites with respect to the population of sites (i.e, all possible locations) 
#' @export
calc_fit_site<-function(sample,pop_df)
{
  
  tot_D = 0
  for(c in c(1:ncol(pop_df)))
  {
    kolmo<-ks.test(sample[,c], ecdf(pop_df[,c]))
    D = kolmo$statistic
    tot_D<-tot_D+D
  }
  
  return(tot_D)
}


#'Function for randomly selecting site subset from all subsets with fitness below a specified percentile. 
#'
#' @param fit_v Vector of fitness scores. 
#' @param prctl Fitness percentile for which any subset of sites with a fitness score below this will be canidates for cross-over.  
#' @return Cumulative kolmogorov smirnov D statistic for each covaraite in the sample of sites with respect to the population of sites (i.e, all possible locations) 
#' @export
select_site<-function(fit_v,prctl)
{
  thresh<-quantile(fit_v,prctl)
  
  canidates<-c(1:length(fit_v))[fit_v<=thresh]
  
  P1<-pop[[sample(canidates,1)]]
  P2<-pop[[sample(canidates,1)]]
  
  parents<-list(P1,P2)
  
  names(parents)<-c('P1','P2')
  
  return(parents)
  
}

#'Function for randomly selecting site subset from all subsets with fitness below a specified percentile. 
#'
#' @param parents List of two parent subsets of sites. 
#' @param mute_rate Rate of mutation 0-1.
#' @param samp_df Data frame containing all covarites of interest as columns, and all potential sites as rows. 
#' @return A new subset of sites that is the random combination of the two provided parent subsets, with mutation. 
#' @export
crossover_site<-function(parents,mute_rate,samp_df)
{
  comb<-c(parents$P1,parents$P2)
  
  offspr<-sample(comb,length(comb)/2)
  
  
  for(i in c(1:length(offspr)))
  {
    rand<-runif(1)
    
    if(rand<mute_rate)
    {
      offspr[i]<-sample(c(1:nrow(samp_df)),1)
    }
    
  }
  
  return(offspr)
}



#'Function for evolving a subset of site locations from a population of canidate site locations for observational regression study design. 
#'
#' @param N Total number of random subsets of sites to initialize,i.e., population size.
#' @param num_sites Maximum number of sites to include in a subset, i.e., total number of sample sites.
#' @param canidate_df Data frame containing all covarites of interest as columns, and all canidate sites locations as rows.
#' @param population_df Data frame containing all covarites of interest as columns, and a discrete approximation of all possible sites locations in the area of interest as rows. 
#' @param generations Total number of generations to iterate. 
#' @param prctl_thresh Fitness percentile for which any subset of sites with a fitness score below this will be canidates for cross-over.
#' @param mute_rate Rate of mutation during cross-over (0-1).
#' @return A new subset of sites that is the random combination of the two provided parent subsets, with mutation. 
#' @export
evolve_site<-function(N,num_sites,canidate_df,population_df,generations,prctl_thresh,mute_rate)
{
  
  pop<-init_pop_site(N,num_sites,canidate_df)
  
  fit_vec<-c()
  
  for(i in c(1:length(pop)))
  {
    fit_vec[i]<-calc_fit_site(target_df[pop[[i]],],target_df)
  }
  
  fit_g<-c()
  
  for(g in c(1:generations))
  {
    
    parents<-select_site(fit_vec,prctl_thresh)
    
    offspring<-crossover_site(parents,mute_rate,canidate_df)
    
    offspring_fit<-calc_fit_site(obs_df[offspring,],pop_df)
    
    min_fit<-which.min(fit_vec)
    max_fit<-which.max(fit_vec)
    
    if(offspring_fit<=fit_vec[min_fit])
    {
      pop[[max_fit]]<-offspring
      fit_vec[max_fit]<-offspring_fit
    }
    
    fit_g[g]<-mean(fit_vec)
    
  }
  
  result_lst <- list(pop,fit_vec,fit_g)
  
  return(result_lst)
}


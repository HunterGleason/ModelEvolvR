


#'Initializes a population of random models.
#'
#' This function initializes a population of random models
#' provided predictors and a response variable as character vectors.
#'
#' @param Y_NAME String, name of response variable in data frame.
#' @param X_NAMES Character vector of possible predictor names present in data frame.
#' @param MAX_TERMS Maximum number of terms to include in a model, does not account for interactions.
#' @param N Total number of random individual models to initialize,i.e., population size.
#' @return A list 'population' of random models as formulas.
#' @export
init_pop<-function(Y_NAME,X_NAMES,MAX_TERMS,N)
{
  pop<-list()
  for(e in c(1:N))
  {
    terms<-sample(c(1:MAX_TERMS),1)

    X_SAMP<-sample(X_NAMES,terms)

    f<-as.formula(paste(Y_NAME,"~",paste(X_SAMP, collapse = '+'),sep=""))

    pop[[e]]<-f
  }

  return(pop)
}


#'Interface to user defined fitness function.
#'
#' This function serves as an interface to user defined fitness function.
#' Provided a formula, or model, and data frame, the function passes
#' calls the users fitness function, which is assumed to only have two
#' parameters, a formula and corresponding data frame,
#' and return a single numeric value 'fitness'.
#'
#' @param FORMULA A formula assumed to correspond to fields present in the data frame.
#' @param DF A data frame containing both the response and predictor variables.
#' @param FIT_FUNC A function that accepts two arguments, a formula and data frame,
#' and returns a fitness as a single numeric value.
#' @return Computed fitness of user defined fitness function.
#' @export
calc_fit<-function(FORMULA, DF, FIT_FUNC)
{
  return(FIT_FUNC(FORMULA,DF))
}

#'Function for selecting two parents from a population of models.
#'
#' This function randomly selects to models 'parents' from a population
#' based on a fitness threshold, i.e., from the pool of models
#' with a fitness score greater than 'or less than' a
#' specified percentile.
#'
#' @param POP The population list.
#' @param FIT_VECT A vector of fitness scores corresponding to the population list.
#' @param PRCTL The percentile at which the fitness threshold is defined for selection.
#' @param MAXIMIZE Boolean, is fitness being maximized, or minimized 'e.g., AIC'
#' @return Two selected parent models in a list.
#' @export
select_parents<-function(POP,FIT_VECT,PRCTL,MAXIMIZE)
{
  parents<-list()

  if(MAXIMIZE==TRUE)
  {
    thresh<-quantile(FIT_VECT,PRCTL, na.rm = T)

    cand_pars<-POP[FIT_VECT>=thresh]

    P1<-sample(cand_pars,1)
    P2<-sample(cand_pars,1)

    parents[[1]]<-P1
    parents[[2]]<-P2
  }else{
    thresh<-quantile(FIT_VECT,PRCTL)

    cand_pars<-POP[FIT_VECT<=thresh]

    P1<-sample(cand_pars,1)
    P2<-sample(cand_pars,1)

    parents[[1]]<-P1
    parents[[2]]<-P2
  }

  return(parents)

}

#'Function for crossing two parent models.
#'
#' This function randomly combines two models, with
#' somewhere between 1 and MAX_TERMS terms. Models are
#' subjected to mutation controlled by a mutation rate.
#'
#' @param PARENTS A list of parent models.
#' @param MUTE_RATE A rate of mutation,0-1.
#' @param MAX_TERMS Maximum number of terms to include in a model, does not account for interactions.
#' @param Y_NAME String name of response variable in data frame.
#' @param X_NAMES Character vector of possible predictors present in data frame.
#' @return A model with genetic material 'input features' from parents, as a formula.
#' @export
cross_over<-function(PARENTS,MUTE_RATE,MAX_TERMS,Y_NAME,X_NAMES)
{

  p1<-all.vars(PARENTS[[1]][[1]])[-1]
  p2<-all.vars(PARENTS[[2]][[1]])[-1]

  comb<-unique(c(p1,p2))

  terms<-sample(c(1:MAX_TERMS),1)

  if(terms>length(comb))
  {

    for(i in c(1:length(comb)))
    {
      rand<-runif(1)

      if(rand<MUTE_RATE)
      {
        comb[i]<-sample(X_NAMES,1)
      }
    }

    return(as.formula(paste(Y_NAME,"~",paste(comb, collapse = '+'),sep="")))
  }else{

    mod<-sample(comb,terms)

    for(i in c(1:length(mod)))
    {
      rand<-runif(1)

      if(rand<MUTE_RATE)
      {
        mod[i]<-sample(X_NAMES,1)
      }
    }

    return(as.formula(paste(Y_NAME,"~",paste(mod, collapse = '+'),sep="")))

  }

}


#'Primary function, performs feature selection through abstraction of evolution process.
#'
#' This function performs feature selection through an abstraction of evolution by
#' minimizing or maximizing a user defined fitness function. The user defined fitness function
#' must take a formula and corresponding data frame as inputs, and return
#' one numeric value representing fitness. The number of generations is up
#' to the user to define, and optimize, however, mean fitness is returned
#' as a vector for each generation to help decide an appropriate number
#' of iterations. This function adapts the mutation rate based on the
#' rate of improvement from one generation to the next in population mean fitness,
#' however a maximum rate of mutation must be provided. This is the rate of mutation
#' when there is no instantaneous change in population mean fitness.
#'
#' @param Y_NAME String name of response variable in data frame.
#' @param X_NAMES Character vector of possible predictors present in data frame.
#' @param DF A data frame containing both the response and predictor variables.
#' @param FIT_FUNC A user defined function that accepts two arguments, a formula and corresponding data frame, and returns a fitness value as a single numeric value.
#' @param MAX_TERMS Maximum number of terms to include in a model, does not account for interactions.
#' @param N Total number of random individual models to initialize,i.e., population size.
#' @param PRCTL The percentile at which the fitness threshold is defined for selection.
#' @param MX_MUTE_RATE A maximum rate of mutation '0-1', i.e., when absolute % slope change in mean fitness equals zero.
#' @param GENZ Total number of generations to be iterated.
#' @param MAXIMIZE Boolean, is fitness being maximized, or minimized 'e.g., AIC'
#' @return A list, with the evolved population list, final fitness vector, and mean fitness vector, in that order.
#' @export
evolve<-function(Y_NAME,X_NAMES,DF,FIT_FUNC,MAX_TERMS,N,PRCTL,MX_MUTE_RATE,GENZ,MAXIMIZE)
{

  set.seed(42)
  pop<-init_pop(Y_NAME,X_NAMES,MAX_TERMS,N)
  fit_vec<-c()

  for(e in c(1:length(pop)))
  {
    fit_vec[e]<-calc_fit(pop[[e]],DF,FIT_FUNC)
  }

  mean_fit<-c(mean(fit_vec))

  MUTE_RATE<-0.0

  for(g in c(1:GENZ))
  {

    parents<-select_parents(pop,fit_vec,PRCTL,MAXIMIZE)

    offspring<-cross_over(parents,MUTE_RATE,MAX_TERMS,Y_NAME,X_NAMES)

    off_fit<-calc_fit(offspring,DF,FIT_FUNC)

    if(MAXIMIZE==TRUE)
    {
      min_fit<-min(fit_vec)
      min_idx<-which.min(fit_vec)

      if(off_fit>min_fit)
      {
        pop[[min_idx]]<-offspring

        fit_vec[min_idx]<-off_fit

      }
    }else{
      max_fit<-max(fit_vec)
      max_idx<-which.max(fit_vec)

      if(off_fit<max_fit)
      {
        pop[[max_idx]]<-offspring

        fit_vec[max_idx]<-off_fit

      }
    }

    mean_fit[g]<-mean(fit_vec)

    if(g>1)
    {
      y1<-mean_fit[g]
      y2<-mean_fit[g-1]
      x1<-g
      x2<-g-1

      slope<-(y1-y2)/(x1-x2)

      p_slope<-abs(slope)*100

      coef_a<-MX_MUTE_RATE/100

      MUTE_RATE = (-coef_a*p_slope) + MX_MUTE_RATE

      if(MUTE_RATE<0)
      {
        MUTE_RATE<-0
      }
    }


  }

  return(list(pop,fit_vec,mean_fit))
}


#Example of fitness function
# fit_fun<-function(FORMULA,DF)
# {
#
#   mod<-summary(lm(formula = FORMULA, data = DF))
#
#   return(mod$adj.r.squared)
#
# }



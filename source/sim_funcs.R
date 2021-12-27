## ---------- Functions to simulate variant detection model ---------- ##

##' Function to set up parameters into desired format
##'
##' @param num_vars number of variants to track (>=2)
##' @param alpha probability an infection of each variant is asymptomatic
##' @param beta_a probability of testing asymptomatic infections
##' @param beta_s probability of testing symptomatic infections
##' @param phi probability that a tested infection of each variant results in a positive test
##' @param gamma probability that a detected infection meets a quality threshold
##' @param omega probability that a sequenced sample produces a high quality genome

setup.params <- function(num_vars,alpha,beta_a,beta_s,phi,gamma,omega){
  
  if(num_vars < 2) {stop("at least two variants must be considered")}
  
  # check that all parameters are either fixed or of length num_vars
  if (length(alpha)!=1 & length(alpha)!=num_vars) 
  {stop("alpha must be numeric or a vector of length num_vars")}
  if (length(phi)!=1 & length(phi)!=num_vars) 
  {stop("phi must be numeric or a vector of length num_vars")}
  if (length(gamma)!=1 & length(gamma)!=num_vars) 
  {stop("gamma must be numeric or a vector of length num_vars")}
  if (length(omega)!=1 & length(omega)!=num_vars) 
  {stop("omega must be numeric or a vector of length num_vars")}
  if (length(beta_a)!=1) {stop("beta_a must be numeric")}
  if (length(beta_s)!=1) {stop("beta_s must be numeric")}
  
  # if a single value is provided, assume this parameter is used for all variants
  if (length(alpha)!=num_vars) {alpha <- rep(alpha,num_vars)}
  if (length(beta_a)!=num_vars) {beta_a <- rep(beta_a,num_vars)}
  if (length(beta_s)!=num_vars) {beta_s <- rep(beta_s,num_vars)}
  if (length(phi)!=num_vars) {phi <- rep(phi,num_vars)}
  if (length(gamma)!=num_vars) {gamma <- rep(gamma,num_vars)}
  if (length(omega)!=num_vars) {omega <- rep(omega,num_vars)}
  
  # put all parameters into a matrix for easy access
  params <- rbind(alpha,beta_a,beta_s,phi,gamma,omega)
  rownames(params) <- c("alpha","beta_a","beta_s","phi","gamma","omega")
  
  return(params)
  
}


##' Function to calculate estimated number of high quality infections given
##' true proportions of variants in the population
##'
##' @param params matrix of relevant parameters
##' @param N infected population size
##' @param P proportion of population infected with each variant

calc.hq.infections <- function(params,N,P){
  
  # confirm that number of variants matches
  if (length(P) != dim(params)[2]) {stop("inconsistent number of variants to evaluate")}
  
  # confirm that proportions add up to one
  if (sum(P) != 1) {stop("variant proportions do not add to one")}
  
  # calculate H: number of detected, high quality samples caused by variant 
  C <- params["phi",] * params["gamma",] *
    (params["alpha",]*params["beta_a",] + ((1-params["alpha",])*params["beta_s",]))
  
  # return estimated number of high quality samples given proportion of population with given variant
  # round to nearest whole number to represent individual samples
  return(round(N*P*C))
  
}


##' Function to calculate estimated number of high quality infections given
##' true proportions of variants in the population
##' with process stochasticity
##'
##' @param params matrix of relevant parameters
##' @param N infected population size
##' @param P proportion of population infected with each variant

calc.hq.infections.rdm <- function(params,N,P){
  
  # confirm that number of variants matches
  if (length(P) != dim(params)[2]) {stop("inconsistent number of variants to evaluate")}
  
  # confirm that proportions add up to one
  if (sum(P) != 1) {stop("variant proportions do not add to one")}
  
  # calculate H: number of detected, high quality samples caused by variant
  
  Nv_a <- rbinom(n=length(P),size=round(N*P),prob=params["alpha",])
  Nv_s <- (round(N*P))-Nv_a
  
  Tv_a <- rbinom(n=length(P),size=Nv_a,prob=params["beta_a",])
  Tv_s <- rbinom(n=length(P),size=Nv_s,prob=params["beta_s",])
  
  Dv_a <- rbinom(n=length(P),size=Tv_a,prob=params["phi",])
  Dv_s <- rbinom(n=length(P),size=Tv_s,prob=params["phi",])
  
  D <- Dv_a + Dv_s
  
  H <- rbinom(n=length(P),size=D,prob=params["gamma",])
  
  # return estimated number of high quality samples given proportion of population with given variant
  # no rounding needed because binom always results in whole numbers
  return(H)
  
}


##' Function to calculate estimated variant numbers given sample size
##'
##' @param params matrix of relevant parameters
##' @param H estimated proportion of high quality samples per variant
##' @param samplesize number of samples that can be sequenced

calc.seq.infections <- function(params,H,samplesize){
  
  # create vector based on H
  samples <- rep(1:length(H),H)
  
  # get total number of samples available to sequence
  tot_samples <- length(samples)
  
  # sample size should not be greater than number of samples available to sequence
  if (samplesize>tot_samples) {
    warning(paste("sample size is larger than estimated number of samples (",
                  as.character(tot_samples),") given parameters",sep=""))
    
    # use total number of available samples for proportion estimation
    samplesize <- tot_samples
  }
  
  # sample from the available samples without replacement
  selected <- sample(samples,size=samplesize,replace=FALSE)
  S <- as.vector(table(factor(selected, levels = 1:length(H)))) # table orders variants by number
  
  # get number of high quality genomes of each variant
  # by multiplying sequenced samples (S) by vector omega
  G <- round(S * params["omega",])
  
  return(G)
  
}


##' Function to calculate estimated variant numbers given sample size
##' with process stochasticity
##'
##' @param params matrix of relevant parameters
##' @param H estimated proportion of high quality samples per variant
##' @param samplesize number of samples that can be sequenced

calc.seq.infections.rdm <- function(params,H,samplesize){
  
  # create vector based on H
  samples <- rep(1:length(H),H)
  
  # get total number of samples available to sequence
  tot_samples <- length(samples)
  
  # sample size should not be greater than number of samples available to sequence
  if (samplesize>tot_samples) {
    warning(paste("sample size is larger than estimated number of samples (",
                  as.character(tot_samples),") given parameters",sep=""))
    
    # use total number of available samples for proportion estimation
    samplesize <- tot_samples
  }
  
  # sample from the available samples without replacement
  selected <- sample(samples,size=samplesize,replace=FALSE)
  S <- as.vector(table(factor(selected, levels = 1:length(H)))) # table orders variants by number
  
  # get number of high quality genomes of each variant
  # by selecting from sequences with probability of omega
  G <- rbinom(n=length(H),size=S,prob=params["omega",])
  
  return(G)
  
}


##' Function to simulate model with two variables
##'
##' @param nsim number of simulations per proportion and sample size
##' @param params matrix of relevant parameters
##' @param N infected population size
##' @param v1 variant proportions to test
##' @param rhos sample proportions to explore

simulate.model.singlevar <- function(nsim,params,N,v1,rhos){
  
  # set up variant proportions as a list of vectors
  v2 <- 1-v1
  Plist <- mapply(c, v1,v2, SIMPLIFY=FALSE)
  
  # get number of available samples
  Hlist <- lapply(Plist,calc.hq.infections,params=params,N=N)
  
  # set up input for fast simulation
  d <- expand.grid(Hlist,rhos)
  colnames(d) <- c("H","rho")
  
  # get sample size from rho and actual H value
  d$samplesize <- apply(d,1,function(x) round(sum(x["H"][[1]])*x["rho"][[1]]))
  
  d <- d[rep(rownames(d),nsim),] # use rep function while doing indexing 
  rownames(d)<-1:nrow(d)
  d$P.in <- rep(Plist,nsim) # add input variant proportions
  
  # calculate output variant proportions
  G <- lapply(1:nrow(d),function(x)
    calc.seq.infections.rdm(params=params,H=d[x,"H"][[1]],samplesize=d[x,"samplesize"]))
  
  outvar <- lapply(G,function(x) x/sum(x))
  
  # add output variant proportions to table for easy access
  d$G <- G
  d$P.out <- outvar
  d$P.diff <- sapply(1:nrow(d),function(x) abs(d[x,"P.in"][[1]][1]-d[x,"P.out"][[1]][1]))
  d$v1.in <- sapply(1:nrow(d),function(x) d[x,"P.in"][[1]][1])
  d$v1.out <- sapply(1:nrow(d),function(x) d[x,"P.out"][[1]][1])
  d$v1.G <- sapply(1:nrow(d),function(x) d[x,"G"][[1]][1])
  
  return(d)
  
}

##' Function to simulate model with two variables
##' with process stochasticity
##'
##' @param nsim number of simulations per proportion and sample size
##' @param params matrix of relevant parameters
##' @param N infected population size
##' @param v1 variant proportions to test
##' @param rhos sample proportions to explore

simulate.model.singlevar.rdm <- function(nsim,params,N,v1,rhos){
  
  # set up variant proportions as a list of vectors
  v2 <- 1-v1
  Plist <- mapply(c, v1,v2, SIMPLIFY=FALSE)
  
  # get number of available samples
  Hlist <- lapply(Plist,calc.hq.infections.rdm,params=params,N=N)
  
  # set up input for fast simulation
  d <- expand.grid(Hlist,rhos)
  colnames(d) <- c("H","rho")
  
  # get sample size from rho and actual H value
  d$samplesize <- apply(d,1,function(x) round(sum(x["H"][[1]])*x["rho"][[1]]))
  
  d <- d[rep(rownames(d),nsim),] # use rep function while doing indexing 
  rownames(d)<-1:nrow(d)
  d$P.in <- rep(Plist,nsim) # add input variant proportions
  
  # calculate output variant proportions
  G <- lapply(1:nrow(d),function(x)
    calc.seq.infections.rdm(params=params,H=d[x,"H"][[1]],samplesize=d[x,"samplesize"]))
  
  outvar <- lapply(G,function(x) x/sum(x))
  
  # add output variant proportions to table for easy access
  d$G <- G
  d$P.out <- outvar
  d$P.diff <- sapply(1:nrow(d),function(x) abs(d[x,"P.in"][[1]][1]-d[x,"P.out"][[1]][1]))
  d$v1.in <- sapply(1:nrow(d),function(x) d[x,"P.in"][[1]][1])
  d$v1.out <- sapply(1:nrow(d),function(x) d[x,"P.out"][[1]][1])
  d$v1.G <- sapply(1:nrow(d),function(x) d[x,"G"][[1]][1])
  
  return(d)
  
}
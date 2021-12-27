## ---------- Sample size calculation functions ---------- ##

# C ratio refers to C1/C2
# Not the inverse C ratio (C2/C1) used in equations
# Functions provided are specific for a two-variant system


# ---------- BIAS FUNCTIONS ---------- #

##' Function to calculate multiplicative bias (observed / actual) in variant prevalence
##'
##' @param p_v1 actual variant prevalence
##' @param c_ratio coefficient of detection ratio

calc.expected.mbias <- function(p_v1,c_ratio){
  bias <- p_v1 + ((1/c_ratio) * (1 - p_v1))
  return(1/bias)
}

##' Function to calculate observed variant prevalence
##'
##' @param p_v1 actual variant prevalence
##' @param c_ratio coefficient of detection ratio

calc.observed.freq <- function(p_v1,c_ratio){
  obs <- p_v1 / (p_v1 + (1/c_ratio)*(1-p_v1))
  return(obs)
}


# ---------- LOGISTIC GROWTH FUNCTIONS ---------- #

##' Function to calculate observed variant prevalence at time t given logistic growth (PDF)
##'
##' @param t time step number (e.g., days) at which to calculate prevalence
##' @param t0 initial variant prevalence (# introductions / population size)
##' @param r logistic growth rate
##' @param c_ratio coefficient of detection ratio, default = 1 (no bias)

calc.freq.logistic <- function(t,t0,r,c_ratio=1){
  a <- (1/t0)-1
  b <- a*(1/c_ratio)
  g <- 1 / (1 + (b*exp(-r*t)))
  return(g)
}

##' Function to cumulative observed variant prevalence at time t given logistic growth (CDF)
##'
##' @param t time step number (e.g., days) at which to calculate CDF
##' @param t0 initial variant prevalence (# introductions / population size)
##' @param r logistic growth rate
##' @param c_ratio coefficient of detection ratio, default = 1 (no bias)

calc.cdf.logistic <- function(t,t0,r,c_ratio=1){
  a <- (1/t0)-1
  b <- a*(1/c_ratio)
  G <- (1/r) * log(b+exp(r*t))
}


# ---------- DETECTION FUNCTIONS ---------- #

##' Function to calculate sample size needed for variant detection
##' assuming cross-sectional sampling
##'
##' @param p_v1 variant prevalence
##' @param prob desired probability of detection
##' @param c_ratio coefficient of detection ratio

calc.samplesize.detect <- function(p_v1,prob,c_ratio){
  p_star <- calc.observed.freq(p_v1,c_ratio)
  n = (log(1-prob)) / (log(1-p_star))
  return(n)
}

##' Function to calculate sample size needed for variant detection
##' assuming periodic sampling
##'
##' @param prob desired probability of detection
##' @param t time step number (e.g., days) at which variant should be detected by
##' @param t0 initial variant prevalence (# introductions / population size)
##' @param r logistic growth rate
##' @param c_ratio coefficient of detection ratio, default = 1 (no bias)

calc.samplesize.detect.cont <- function(prob,t,t0,r,c_ratio=1){
  n = -(log(1-prob)) / (calc.cdf.logistic(t,t0,r,c_ratio)-calc.cdf.logistic(0,t0,r,c_ratio))
  return(n)
}

##' Function to probability of detecting a variant given a per-timestep sample size
##' assuming periodic sampling
##'
##' @param n per-timestep (e.g., per day) sample size
##' @param t time step number (e.g., days) at which variant should be detected by
##' @param t0 initial variant prevalence (# introductions / population size)
##' @param r logistic growth rate
##' @param c_ratio coefficient of detection ratio, default = 1 (no bias)

calc.prob.detect.cont <- function(n,t,t0,r,c_ratio=1){
  prob <- 1-exp(-n*(calc.cdf.logistic(t,t0,r,c_ratio)-calc.cdf.logistic(0,t0,r,c_ratio)))
  return(prob)
}


# ---------- PREVALENCE MEASURING FUNCTIONS ---------- #

##' Function to calculate sample size needed for variant prevalence estimation
##' assuming cross-sectional sampling
##'
##' @param p_v1 variant prevalence
##' @param prob desired confidence in variant prevalence estimate
##' @param precision desired precision in variant prevalence estimate
##' @param c_ratio coefficient of detection ratio

calc.samplesize.prev <- function(p_v1,prob,precision,c_ratio){
  p_star <- calc.observed.freq(p_v1,c_ratio)
  Z <- qnorm(1-((1-prob)/2))
  n = ((Z^2)*p_star*(1-p_star)) / ((p_star*precision)^2)
  return(n)
}
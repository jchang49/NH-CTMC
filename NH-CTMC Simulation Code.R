rm(list = ls())
set.seed(0)

# Load packages for simulation
library('bbmle')
library('zoo')
library('msm')
library('gdata')
library('cubature')

## Set Simulation Parameters
# Number of subjects
n = 100
# Number of detection timepoint (Length of the Markov Chain)
T = 50
# Parameter for lambda (rate of trasitioning from state 0 to state 1)
a = 9
# Parameter for mu (rate of transitioning from state 1 to state 0)
b = 5
# Binary covariate coefficient for transitioning from state 0 to state 1
Beta01 = 1
# Binary covariate coefficient for transitioning from state 1 to state 0
Beta10 = 0.5
# Probability of starting from state 1
p = 0.5
# Probability of being in group 1 of binary covariate
p.covariate = 0.5


## Define Transition Probability Functions
# Probability of remaining at state 0
P_00 = function(a, b, x, Beta01, Beta10, s, t) {
  (exp(x * Beta10) * cubintegrate(function(u){(1 + (u - 1)^(a / b))^exp(x * Beta01) * u^(exp(x * Beta10) - 1)}, 1 + s^b, 1 + (s + t)^b)$integral + (1 + s^a)^(exp(x * Beta01)) * (1 + s^b)^(exp(x * Beta10))) / ((1 + (s + t)^a)^exp(x * Beta01) * (1 + (s + t)^b)^exp(x * Beta10))
}
# Probability of transitioning from state 0 to state 1
P_01 = function(a, b, x, Beta01, Beta10, s, t) {
  1 - P_00(a, b, x, Beta01, Beta10, s, t)
}
# Probability of remaining at state 1
P_11 = function(a, b, x, Beta01, Beta10, s, t) {
  (exp(x * Beta01) * cubintegrate(function(v){(1 + (v - 1)^(b / a))^exp(x * Beta10) * v^(exp(x * Beta01) - 1)}, 1 + s^a, 1 + (s + t)^a)$integral + (1 + s^a)^(exp(x * Beta01)) * (1 + s^b)^(exp(x * Beta10))) / ((1 + (s + t)^a)^exp(x * Beta01) * (1 + (s + t)^b)^exp(x * Beta10))
}
# Probability of transitioning from state 1 to state 0
P_10 = function(a, b, x, Beta01, Beta10, s, t) {
  1 - P_11(a, b, x, Beta01, Beta10, s, t)
}



## Define Transition Rate Functions
# Transition Rate Functions for transitioning from state 0 to state 1
delta1 = function (t) {
  (a * t^(a - 1)) / (1 + t^a) * exp(x*Beta01)
}
# Transition Rate Functions for transitioning from state 1 to state 0
delta2 = function (t) {
  (b * t^(b - 1)) / (1 + t^b) * exp(x*Beta10)
}


## Non-homogeneous continuous time Markov Chain simulation
# Set up empty matrix for storing states of the markov chains
MC = matrix(rep(NA, n * (T + 1)), nrow = n, ncol = (T + 1))
# Generate values for binary covariates for each subject through Bernoulli distribution
covariate = rbinom(n, size = 1, p.covariate)

# Generate non-homogeneous continuous time Markov chain for each subject by simulating non-homogeneous Poisson process (with Thinning process)
# See Algorithm 2
for (i in 1:n) {
  # Set initial state for each subject through Bernoulli distribution
  state = c(rbinom(n = 1,size = 1, p))
  timer = c(0)
  x = covariate[i]
  
  while (timer[length(timer)] <= T) {
    tm = timer[length(timer)]
    # Maximum rate within the time interval (for thinning process)
    DELTA_MAX = (optimize(delta1, interval = c(tm, T), maximum = TRUE)$objective) * as.numeric(state[length(state)] == 0) + 
                (optimize(delta2, interval = c(tm, T), maximum = TRUE)$objective) * as.numeric(state[length(state)] == 1)
    
    # Simulate like homogeneous Poisson process
    # Simulate the number of events with Poisson distribution
    ntemp = rpois(1, lambda = DELTA_MAX * (T - tm))
    # Uniformly distribute the event times throughout time interval. 
    s = runif(ntemp, tm, T)
    s = sort(s)
    p_temp = (delta1(s) * as.numeric(state[length(state)] == 0) + delta2(s) * as.numeric(state[length(state)] == 1)) / DELTA_MAX
    
    # Drop the events through thinning process
    s0 = rbinom(ntemp, size = 1, prob = replace(p_temp, p_temp > 1, 1))
    v = s[s0 == 1]
    
    # If there is no event remaining, there is no transition until the end time point. 
    if (length(v) == 0) break
    
    timer = c(timer, v[1])
    state = c(state, as.numeric(state[length(state)] == 0))
  }
  
  # Record the observations in discrete integer time points. 
  for (j in 1:length(timer)) {
    if (timer[j] > T) {
      break
    }
    else {
      MC[i, ceiling(timer[j]) + 1] = state[j]
    }
  }
}

# Fill in the "NA" blanks with previously observed state
for (i in 1:nrow(MC)) {
  for (j in 1:ncol(MC)) {
    if (is.na(MC[i, j])) {
      MC[i, j] = MC[i, j-1]
    }
  }
}


## Non-homogeneous continuous time Markov Chain Estimation
# Calculate likelihood function of all the Markov chains in the matrix 
nLL = function(a, b, Beta01, Beta10) {
  sum = 0
  for (i in 1:nrow(MC)) {
    for (j in 2:(ncol(MC))) {
      sum = sum - 
        log(P_00(a, b, covariate[i], Beta01, Beta10, j-2, 1)) * (MC[i, j-1]==0) * (MC[i, j]==0) - 
        log(P_01(a, b, covariate[i], Beta01, Beta10, j-2, 1)) * (MC[i, j-1]==0) * (MC[i, j]==1) - 
        log(P_10(a, b, covariate[i], Beta01, Beta10, j-2, 1)) * (MC[i, j-1]==1) * (MC[i, j]==0) - 
        log(P_11(a, b, covariate[i], Beta01, Beta10, j-2, 1)) * (MC[i, j-1]==1) * (MC[i, j]==1)
    }
  }
  return(sum)
}

# Starting values for MLE's optimization process
a.mle.initial = 1
b.mle.initial = 1
Beta.mle.initial = 0

# MLE Calculation using Nelder-Mead method
fit = mle2(minuslogl = nLL, 
           start = list(a = a.mle.initial, b = b.mle.initial, Beta01 = Beta.mle.initial, Beta10 = Beta.mle.initial), 
           optimizer = 'constrOptim', 
           ui = rbind(c(1,0,0,0), 
                      c(0,1,0,0), 
                      c(0,0,0,0), 
                      c(0,0,0,0)), 
           ci = c(0,0,-1,-1), 
           method = 'Nelder-Mead')
fit2 = summary(fit)
fit2

## Utility Functions
# Laste edited: 101722

## (1) buchanan_log10N, (2) gompertz_log10N, (3) baranyi_log10N
## (4) buchanan_without_lag_log10N, (5) baranyi_without_lag_log10N
# Source: All equations for the 5 models below are copied from the Nlms package in R (https://rdrr.io/cran/nlsMicrobio/src/R/growthmodels.R)

# Purpose: Calculate log10N using respective growth model (either buchanan, gompertz, baranyi, buchanan_without_lag, or baranyi_without_lag)

# Parameters: Same for (1) buchanan, (2) gompertz, and (3) baranyi
# (i) t: time in hours
# (ii) lag: length of lag phase
# (iii) mumax: growth rate (ln/hour)
# (iv) LOG10N0: initial microbial concentration 
# (v) LOG10Nmax: carrying capacity

# Parameters: Same for (4) buchanan_without_lag, and (5) baranyi_without_lag
# (i) t: time in hours
# (ii) mumax: growth rate (ln/hour)
# (iii) LOG10N0: initial microbial concentration 
# (iv) LOG10Nmax: carrying capacity

# Functions: 
buchanan_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax) {
  ans <- LOG10N0 + (t >= lag) * (t <= (lag + (LOG10Nmax - LOG10N0) * log(10)/mumax)) * mumax * (t - lag)/log(10) + (t >= lag) * (t > (lag + (LOG10Nmax - LOG10N0) * log(10)/mumax)) * (LOG10Nmax - LOG10N0)
  return(ans)
}

gompertz_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax) {
  ans <- LOG10N0 + (LOG10Nmax - LOG10N0) * exp(-exp(mumax * exp(1) * (lag - t)/((LOG10Nmax - LOG10N0) * log(10)) + 1))
  return(ans)
}

baranyi_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax) {
  ans <- LOG10Nmax + log10((-1 + exp(mumax * lag) + exp(mumax * t))/(exp(mumax * t) - 1 + exp(mumax * lag) * 10^(LOG10Nmax - LOG10N0)))
  return(ans)
}

buchanan_without_lag_log10N = function(t,mumax,LOG10N0,LOG10Nmax) {
  ans <- LOG10N0 + (t <= ((LOG10Nmax - LOG10N0) * log(10)/mumax)) * mumax * t/log(10) + (t > ((LOG10Nmax - LOG10N0) * log(10)/mumax)) * (LOG10Nmax - LOG10N0)
  return(ans)
}

baranyi_without_lag_log10N = function(t,mumax,LOG10N0,LOG10Nmax) {
  ans <- (LOG10Nmax - log10(1 + (10^(LOG10Nmax - LOG10N0) - 1) * exp(-mumax * t)))
  return(ans)
}

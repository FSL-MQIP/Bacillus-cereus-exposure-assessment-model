# Load packages 
library(nlsMicrobio)
library(minpack.lm)

# Import dataset 
data <-read.csv("data_16_new.csv")
colnames(data) = c("t","Rep","LOG10N","Isolate")

# Subset data 
Iso135Rep2 = subset(data,Rep == "rep2" & Isolate == 135)

# Plot 
plot(Iso135Rep2$t,Iso135Rep2$LOG10N)
lines(Iso135Rep2$t,fitted(Iso135Rep2.bar_LM),col="red")
mod<-lm(Iso135Rep2$LOG10N[6:8]~Iso135Rep2$t[6:8])
slope = coef(mod)[2]*log(10)

# Baranyi
baranyi_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax) {
  ans <- LOG10Nmax + log10((-1 + exp(mumax * lag) + exp(mumax * t))/(exp(mumax * t) - 1 + exp(mumax * lag) * 10^(LOG10Nmax - LOG10N0)))
  return(ans)
}

# Fitting Baranyi model 
Iso135Rep2.bar_LM <- nlsLM(LOG10N ~ baranyi_log10N(t,lag,mumax,LOG10N0,LOG10Nmax), data=Iso135Rep2,
                           start=list (
                           LOG10N0 = 3.8,
                           lag = 8, 
                           mumax = 0.389, 
                           LOG10Nmax = 6.3), 
                           lower = c(0,0,0,0))


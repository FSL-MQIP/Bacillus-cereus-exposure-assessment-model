# Load packages 
library(nlsMicrobio)
library(minpack.lm)


# Import data 
# 22dC data 
data1 <-read.csv("Primary model/InputFiles/data_22_new.csv")
# 10dC data 
data2 <- read.csv("Primary model/InputFiles/data_10_new.csv")

# Load primary growth models
source("Primary model/UtilityFunctions_baranyi.R")

# Subset data 
# Iso433Rep1 @ 22dC 
Iso433Rep1 <- subset(data1, Isolate == "433" & Rep =="rep1")
# Iso433Rep2 @ 22dC 
Iso433Rep2 <- subset(data1, Isolate == "433" & Rep =="rep2")
# Iso413Rep1 @ 10dC 
Iso413Rep1 <- subset(data2, Isolate == "413" & Rep =="rep1")
Iso413Rep1_new <- Iso413Rep1[1:6,]
colnames(Iso413Rep1_new) <- c("t","Rep","LOG10N","Isolate")
Iso413Rep1_new$t <- as.numeric(Iso413Rep1_new$t)
Iso413Rep1_new$LOG10N <- as.numeric(Iso413Rep1_new$LOG10N)

# Fit Baranyi 
Iso413Rep1.bar_LM <- nlsLM(LOG10N ~ baranyi_log10N(t,lag,mumax,LOG10N0,LOG10Nmax), data=Iso413Rep1_new,
                           start=list (
                            LOG10N0 = 3.5,
                            lag = 24, 
                            mumax =  0.04986383, 
                            LOG10Nmax = 7),
                            lower = c(0,0,0,0))
                           

# Fitting Baranyi without lag 
Iso433Rep1.bar_nolag_LM<- nlsLM(LOG10N ~ baranyi_without_lag_log10N(t,mumax,LOG10N0,LOG10Nmax), data=Iso433Rep1,
                                start=list (
                                LOG10N0 = 3.5,
                                mumax = 0.3227814, 
                                LOG10Nmax = 7), 
                                lower = c(0,0,0))

Iso433Rep2.bar_nolag_LM<- nlsLM(LOG10N ~ baranyi_without_lag_log10N(t,mumax,LOG10N0,LOG10Nmax), data=Iso433Rep2,
                                start=list (
                                LOG10N0 = 3.5,
                                mumax = 0.3227814, 
                                LOG10Nmax = 7), 
                                lower = c(0,0,0))

# Plot fitted Iso413Rep1 @ 10dC 
plot(Iso413Rep1_new$t,Iso413Rep1_new$LOG10N,
     xlim= c(0,300),
     ylim= c(3,8),
     xlab="t (h)",
     ylab= "LOG10N",
     main="Iso413Rep1_new")

lines(Iso413Rep1_new$t,fitted(Iso413Rep1.bar_LM),col="red")

# Load packages 
library(nlsMicrobio)
library(minpack.lm)
library(stats)


# Import data 
data1 <-read.csv("data1.csv")
data1 <- data1[,2:5]
data2 <-read.csv("data2.csv") #data added for Isolate457Rep3 & Rep4, Isolate536Rep3
data_10dC <- read.csv("SMB_10degC_Compiled.csv")
data_10dC <- data_10dC[,1:4]
colnames(data_10dC) <- c('t','Rep','LOG10N','Isolate')
data_10dC<-na.omit(data_10dC) 


# Load primary growth models
source("UtilityFunctions.R")


#Note: This is an example of Isolate457Rep2. The code applies to other samples by changing sample names

# Subset data 
Iso649Rep2 <- subset(data_10dC, Isolate == "649" & Rep =="rep2")


# Find starting values
plot(Iso649Rep2$t,Iso649Rep2$LOG10N)
mod<-lm(Iso649Rep2$LOG10N[4:7] ~ Iso649Rep2$t[4:7])
slope<-coef(mod)[2]
slope

# Plot fitted growth curve
lines(Iso413Rep1$t,fitted(Iso413Rep1.buc_LM),col="red")
lines(Iso413Rep2$t,fitted(Iso413Rep2.gom_LM),col="blue")
lines(Iso413Rep2$t,fitted(Iso413Rep2.bar_LM),col="green")


# Fit primary growth models 
# Fit Buchanan 
Iso649Rep2.buc_LM<- nlsLM(LOG10N ~ buchanan_log10N(t,lag,mumax,LOG10N0,LOG10Nmax), data=Iso649Rep2,
                          start=list (
                          LOG10N0 = 3.5,
                          lag = 48, 
                          mumax = 0.01467, 
                          LOG10Nmax = 6.5), 
                          lower = c(0,0,0,0))

# Fit Gompertz 
Iso649Rep2.gom_LM <- nlsLM(LOG10N ~ gompertz_log10N(t,lag,mumax,LOG10N0,LOG10Nmax), data=Iso649Rep2,
                           start=list (
                           LOG10N0 = 3.5,
                           lag = 48,
                           mumax = 0.01467,
                           LOG10Nmax = 6.5), 
                           lower = c(0,0,0,0))

# Fit Baranyi 
Iso649Rep2.bar_LM <- nlsLM(LOG10N ~ baranyi_log10N(t,lag,mumax,LOG10N0,LOG10Nmax), data=Iso649Rep2,
                           start=list (
                           LOG10N0 = 3.5,
                           lag = 48, 
                           mumax = 0.01467, 
                           LOG10Nmax = 6.5), 
                           lower = c(0,0,0,0))

# Collect growth parameters estimated from 3 primary growth models
Iso649Rep2_buccoef<-coef(Iso649Rep2.buc_LM)
Iso649Rep2_gomcoef<-coef(Iso649Rep2.gom_LM)
Iso649Rep2_barcoef<-coef(Iso649Rep2.bar_LM)

# Compute AIC, BIC, RMSE for 3 primary growth models
# For Isolate457, Rep3 & 4 are used for RMSE calculation 
# For Isolate536, Rep2 & 3 are used for RMSE calculation  
AIC_Iso649Rep2_buc=AIC(Iso649Rep2.buc_LM)
BIC_Iso649Rep2_buc=BIC(Iso649Rep2.buc_LM)
Iso649resid.buc = c(Iso649Rep1.buc_LM$m$resid(), Iso649Rep2.buc_LM$m$resid())
Iso649rmse.buc = sqrt(sum(Iso649resid.buc^2)/length(Iso649resid.buc))

AIC_Iso649Rep2_gom=AIC(Iso649Rep2.gom_LM)
BIC_Iso649Rep2_gom=BIC(Iso649Rep2.gom_LM)
Iso649resid.gom = c(Iso649Rep1.gom_LM$m$resid(), Iso649Rep2.gom_LM$m$resid())
Iso649rmse.gom = sqrt(sum(Iso649resid.gom^2)/length(Iso649resid.gom))


AIC_Iso649Rep2_bar=AIC(Iso649Rep2.bar_LM)
BIC_Iso649Rep2_bar=BIC(Iso649Rep2.bar_LM)
Iso649resid.bar = c(Iso649Rep1.bar_LM$m$resid(), Iso649Rep2.bar_LM$m$resid())
Iso649rmse.bar = sqrt(sum(Iso649resid.bar^2)/length(Iso649resid.bar))

#Note: Isolate433Rep2 has an estimated lag time of 0 after fitting into Buchanan, Gompertz, Baranyi models
#Isolate433Rep1 & Rep2 are then fitted into Buchanan_without-lag and Baranyi_without_lag

# Fit Buchanan without lag 
Iso193Rep2.buc_nolag_LM<- nlsLM(LOG10N ~ buchanan_without_lag_log10N(t,mumax,LOG10N0,LOG10Nmax), data=Iso193Rep2,
                                start=list (
                                LOG10N0 = 3,
                                mumax = 0.0155, 
                                LOG10Nmax = 5.2), 
                                lower = c(0,0,0))

# Fitting Baranyi without lag 
Iso193Rep2.bar_nolag_LM<- nlsLM(LOG10N ~ baranyi_without_lag_log10N(t,mumax,LOG10N0,LOG10Nmax), data=Iso193Rep2,
                                start=list (
                                LOG10N0 = 3,
                                mumax = 0.0155, 
                                LOG10Nmax = 5.2), 
                                lower = c(0,0,0))

# Collect growth parameters estimated from 2 primary growth models (without lag)
Iso193Rep2_buc_nolag_coef<-coef(Iso193Rep2.buc_nolag_LM)
Iso193Rep2_bar_nolag_coef<-coef(Iso193Rep2.bar_nolag_LM)

# Compute AIC, BIC, RMSE for 2 primary growth models (without lag)
AIC_Iso193Rep2_buc_nolag=AIC(Iso193Rep2.buc_nolag_LM)
BIC_Iso193Rep2_buc_nolag=BIC(Iso193Rep2.buc_nolag_LM)
Iso193resid.buc_nolag = c(Iso193Rep1.buc_nolag_LM$m$resid(), Iso193Rep2.buc_nolag_LM$m$resid())
Iso193rmse.buc_nolag = sqrt(sum(Iso193resid.buc_nolag^2)/length(Iso193resid.buc_nolag))

AIC_Iso193Rep2_bar_nolag=AIC(Iso193Rep2.bar_nolag_LM)
BIC_Iso193Rep2_bar_nolag=BIC(Iso193Rep2.bar_nolag_LM)
Iso193resid.bar_nolag = c(Iso193Rep1.bar_nolag_LM$m$resid(), Iso193Rep2.bar_nolag_LM$m$resid())
Iso193rmse.bar_nolag = sqrt(sum(Iso193resid.bar_nolag^2)/length(Iso193resid.bar_nolag))

# Combine growth parameters, AIC, BIC data for each isolate 
Iso457parameters<- rbind(Iso457Rep1_buccoef,Iso457Rep1_gomcoef,Iso457Rep1_barcoef,
                         Iso457Rep2_buccoef,Iso457Rep2_gomcoef,Iso457Rep2_barcoef)
AIC<-rbind(AIC_Iso457Rep1_buc,AIC_Iso457Rep1_gom,AIC_Iso457Rep1_bar,
           AIC_Iso457Rep2_buc,AIC_Iso457Rep2_gom,AIC_Iso457Rep2_bar)
rownames(AIC)<-c("Iso457Rep1_buccoef","Iso457Rep1_gomcoef","Iso457Rep1_barcoef",
                 "Iso457Rep2_buccoef","Iso457Rep2_gomcoef","Iso457Rep2_barcoef")
colnames(AIC)<-c("AIC")
BIC<-rbind(BIC_Iso457Rep1_buc,BIC_Iso457Rep1_gom,BIC_Iso457Rep1_bar,
           BIC_Iso457Rep2_buc,BIC_Iso457Rep2_gom,BIC_Iso457Rep2_bar)
rownames(BIC)<-c("Iso457Rep1_buccoef","Iso457Rep1_gomcoef","Iso457Rep1_barcoef",
                 "Iso457Rep2_buccoef","Iso457Rep2_gomcoef","Iso457Rep2_barcoef")
colnames(BIC)<-c("BIC")
Result_457<-cbind(Iso457parameters,AIC,BIC)

# Combine RMSE for each isolate
RMSE_457<-rbind(Iso457rmse.buc,Iso457rmse.gom,Iso457rmse.bar)
colnames(RMSE_457)<-c("RMSE")

#Generate output
Result<-rbind(Result_193,Result_194,Result_125,Result_135,Result_402,Result_407,Result_413,Result_433,
              Result_433_nolag,Result_457,Result_457Rep3,Result_457Rep4,Result_474,Result_495,Result_518,
              Result_536,Result_536Rep3,Result_564,Result_570,Result_638,Result_649)
as.data.frame(Result)
write.csv(Result,"Result.csv")

RMSE<-rbind(RMSE_193,RMSE_194,RMSE_125,RMSE_135,RMSE_402,RMSE_407,RMSE_413,RMSE_433,RMSE_457,
            RMSE_474,RMSE_495,RMSE_518,RMSE_536,RMSE_564,RMSE_570,RMSE_638,RMSE_649)
as.data.frame(RMSE)
write.csv(RMSE,"RMSE.csv")




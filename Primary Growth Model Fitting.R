# Load packages 
library(nlsMicrobio)
library(minpack.lm)
library(stats)

rm(Iso413Rep2.buc_nolag_LM)

# Import data 
data1 <-read.csv("data1.csv")
data1 <- data1[,2:5]
data2 <-read.csv("data2.csv") #data added for Isolate457Rep3 & Rep4, Isolate536Rep3


# Load primary growth models
source("UtilityFunctions.R")


#Note: This is an example of Isolate457Rep2. The code applies to other samples by changing sample names

# Subset data 
Iso457Rep2 <- subset(data1, Isolate == "457" & Rep =="rep2")


# Fit primary growth models 
# Fit Buchanan 
Iso457Rep2.buc_LM<- nlsLM(LOG10N ~ buchanan_log10N(t,lag,mumax,LOG10N0,LOG10Nmax), data=Iso457Rep2,
                          start=list (
                          LOG10N0 = 3.5,
                          lag = 15, 
                          mumax = 0.2945, 
                          LOG10Nmax = 6.3), 
                          lower = c(0,0,0,0))

# Fit Gompertz 
Iso457Rep2.gom_LM <- nlsLM(LOG10N ~ gompertz_log10N(t,lag,mumax,LOG10N0,LOG10Nmax), data=Iso457Rep2,
                           start=list (
                           LOG10N0 = 4.2,
                           lag = 18,
                           mumax = 0.3272,
                           LOG10Nmax = 6.227), 
                           lower = c(0,0,0,0))

# Fit Baranyi 
Iso457Rep2.bar_LM <- nlsLM(LOG10N ~ baranyi_log10N(t,lag,mumax,LOG10N0,LOG10Nmax), data=Iso457Rep2,
                           start=list (
                           LOG10N0 = 3.5,
                           lag = 14, 
                           mumax = 0.4, 
                           LOG10Nmax = 6), 
                           lower = c(0,0,0,0))

# Collect growth parameters estimated from 3 primary growth models
Iso457Rep2_buccoef<-coef(Iso457Rep2.buc_LM)
Iso457Rep2_gomcoef<-coef(Iso457Rep2.gom_LM)
Iso457Rep2_barcoef<-coef(Iso457Rep2.bar_LM)

# Compute AIC, BIC, RMSE for 3 primary growth models
# For Isolate457, Rep3 & 4 are used for RMSE calculation 
# For Isolate536, Rep2 & 3 are used for RMSE calculation  
AIC_Iso457Rep2_buc=AIC(Iso457Rep2.buc_LM)
BIC_Iso457Rep2_buc=BIC(Iso457Rep2.buc_LM)
Iso457resid.buc = c(Iso457Rep3.buc_LM$m$resid(), Iso457Rep4.buc_LM$m$resid())
Iso457rmse.buc = sqrt(sum(Iso457resid.buc^2)/length(Iso457resid.buc))

AIC_Iso457Rep2_gom=AIC(Iso457Rep2.gom_LM)
BIC_Iso457Rep2_gom=BIC(Iso457Rep2.gom_LM)
Iso457resid.gom = c(Iso457Rep3.gom_LM$m$resid(), Iso457Rep4.gom_LM$m$resid())
Iso457rmse.gom = sqrt(sum(Iso457resid.gom^2)/length(Iso457resid.gom))


AIC_Iso457Rep2_bar=AIC(Iso457Rep2.bar_LM)
BIC_Iso457Rep2_bar=BIC(Iso457Rep2.bar_LM)
Iso457resid.bar = c(Iso457Rep3.bar_LM$m$resid(), Iso457Rep4.bar_LM$m$resid())
Iso457rmse.bar = sqrt(sum(Iso457resid.bar^2)/length(Iso457resid.bar))

#Note: Isolate433Rep2 has an estimated lag time of 0 after fitting into Buchanan, Gompertz, Baranyi models
#Isolate433Rep1 & Rep2 are then fitted into Buchanan_without-lag and Baranyi_without_lag

# Fit Buchanan without lag 
Iso433Rep2.buc_nolag_LM<- nlsLM(LOG10N ~ buchanan_without_lag_log10N(t,mumax,LOG10N0,LOG10Nmax), data=Iso433Rep2,
                                start=list (
                                LOG10N0 = 3.2,
                                mumax = 0.058, 
                                LOG10Nmax = 6.2), 
                                lower = c(0,0,0))

# Fitting Baranyi without lag 
Iso433Rep2.bar_nolag_LM<- nlsLM(LOG10N ~ baranyi_without_lag_log10N(t,mumax,LOG10N0,LOG10Nmax), data=Iso433Rep2,
                                start=list (
                                LOG10N0 = 3.5,
                                mumax = 0.158, 
                                LOG10Nmax = 7), 
                                lower = c(0,0,0))

# Collect growth parameters estimated from 2 primary growth models (without lag)
Iso433Rep2_buc_nolag_coef<-coef(Iso433Rep2.buc_nolag_LM)
Iso433Rep2_bar_nolag_coef<-coef(Iso433Rep2.bar_nolag_LM)

# Compute AIC, BIC, RMSE for 2 primary growth models (without lag)
AIC_Iso433Rep2_buc_nolag=AIC(Iso433Rep2.buc_nolag_LM)
BIC_Iso433Rep2_buc_nolag=BIC(Iso433Rep2.buc_nolag_LM)
Iso433resid.buc_nolag = c(Iso433Rep1.buc_nolag_LM$m$resid(), Iso433Rep2.buc_nolag_LM$m$resid())
Iso433rmse.buc_nolag = sqrt(sum(Iso433resid.buc_nolag^2)/length(Iso433resid.buc_nolag))

AIC_Iso433Rep2_bar_nolag=AIC(Iso433Rep2.bar_nolag_LM)
BIC_Iso433Rep2_bar_nolag=BIC(Iso433Rep2.bar_nolag_LM)
Iso433resid.bar_nolag = c(Iso433Rep1.bar_nolag_LM$m$resid(), Iso433Rep2.bar_nolag_LM$m$resid())
Iso433rmse.bar_nolag = sqrt(sum(Iso433resid.bar_nolag^2)/length(Iso433resid.bar_nolag))

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




# Load packages 
library(nlsMicrobio)
library(minpack.lm)
library(stats)


# Import data 
# 22dC data 
data1 <-read.csv("data1.csv")
data1 <- data1[,2:5]
data2 <-read.csv("data2.csv") #data added for Isolate457Rep3 & Rep4, Isolate536Rep3

# 10dC data 
data_10dC <- read.csv("SMB_10degC_Compiled.csv")
data_10dC <- data_10dC[,1:4]
data_10dC_new <- read.csv("Iso_193_10degC.csv")


# Load primary growth models
source("UtilityFunctions.R")


#Note: This is an example of Isolate193Rep2 (10dC data). The code applies to other samples by changing sample names

# Subset data 
Iso413Rep1_10dC <- subset(data_10dC, Isolate == "413" & Rep =="rep1")
colnames(Iso413Rep1_10dC) = c("t","Rep","LOG10N","Isolate")

plot(Iso413Rep1_10dC$t,Iso413Rep1_10dC$LOG10N)
mod <- lm(Iso649Rep2$LOG10N ~ Iso649Rep2$t)
slope <- coef(mod) [2]
slope


# Fit primary growth models 
# Fit Buchanan 
Iso193Rep2.buc_LM<- nlsLM(LOG10N ~ buchanan_log10N(t,lag,mumax,LOG10N0,LOG10Nmax), data=Iso193Rep2,
                          start=list (
                          LOG10N0 = 3,
                          lag = 100, 
                          mumax = 0.006, 
                          LOG10Nmax = 6), 
                          lower = c(0,0,0,0))

# Fit Gompertz 
Iso193Rep2.gom_LM <- nlsLM(LOG10N ~ gompertz_log10N(t,lag,mumax,LOG10N0,LOG10Nmax), data=Iso193Rep2,
                           start=list (
                           LOG10N0 = 3,
                           lag = 100,
                           mumax = 0.006,
                           LOG10Nmax = 6), 
                           lower = c(0,0,0,0))

# Fit Baranyi 
Iso649Rep2.bar_LM <- nlsLM(LOG10N ~ baranyi_log10N(t,lag,mumax,LOG10N0,LOG10Nmax), data=Iso649Rep2,
                           start=list (
                           LOG10N0 = 3,
                           lag = 8, 
                           mumax = 0.43, 
                           LOG10Nmax = 6.5), 
                           lower = c(0,0,0,0))

# Collect growth parameters estimated from 3 primary growth models
Iso193Rep2_buccoef<-coef(Iso193Rep2.buc_LM)
Iso193Rep2_gomcoef<-coef(Iso193Rep2.gom_LM)
Iso193Rep2_barcoef<-coef(Iso193Rep2.bar_LM)

# Compute AIC, BIC, RMSE for 3 primary growth models
AIC_Iso193Rep2_buc=AIC(Iso193Rep2.buc_LM)
BIC_Iso193Rep2_buc=BIC(Iso193Rep2.buc_LM)
Iso193resid.buc = c(Iso193Rep1.buc_LM$m$resid(), Iso193Rep2.buc_LM$m$resid())
Iso193rmse.buc = sqrt(sum(Iso193resid.buc^2)/length(Iso193resid.buc))

AIC_Iso193Rep1_gom=AIC(Iso193Rep1.gom_LM)
BIC_Iso193Rep1_gom=BIC(Iso193Rep1.gom_LM)
Iso193resid.gom = c(Iso193Rep1.gom_LM$m$resid(), Iso193Rep2.gom_LM$m$resid())
Iso193rmse.gom = sqrt(sum(Iso193resid.gom^2)/length(Iso193resid.gom))


AIC_Iso193Rep1_bar=AIC(Iso193Rep1.bar_LM)
BIC_Iso193Rep1_bar=BIC(Iso193Rep1.bar_LM)
Iso193resid.bar = c(Iso193Rep1.bar_LM$m$resid(), Iso193Rep2.bar_LM$m$resid())
Iso193rmse.bar = sqrt(sum(Iso193resid.bar^2)/length(Iso193resid.bar))

#Note: Isolate433Rep2 (22dC data) has an estimated lag time of 0 after fitting into Buchanan, Gompertz, Baranyi models
#Isolate433Rep1&2 are then fitted into Buchanan_without-lag and Baranyi_without_lag

# Fit Buchanan without lag 
Iso433Rep2.buc_nolag_LM<- nlsLM(LOG10N ~ buchanan_without_lag_log10N(t,mumax,LOG10N0,LOG10Nmax), data=Iso433Rep2,
                                start=list (
                                LOG10N0 = 3.5,
                                mumax = 0.3, 
                                LOG10Nmax = 7), 
                                lower = c(0,0,0))

# Fitting Baranyi without lag 
Iso433Rep2.bar_nolag_LM<- nlsLM(LOG10N ~ baranyi_without_lag_log10N(t,mumax,LOG10N0,LOG10Nmax), data=Iso433Rep2,
                                start=list (
                                LOG10N0 = 3.5,
                                mumax = 0.3, 
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
Iso193parameters<- rbind(Iso193Rep1_buccoef,Iso193Rep1_gomcoef,Iso193Rep1_barcoef,
                         Iso193Rep2_buccoef,Iso193Rep2_gomcoef,Iso193Rep2_barcoef)
AIC<-rbind(AIC_Iso193Rep1_buc,AIC_Iso193Rep1_gom,AIC_Iso193Rep1_bar,
           AIC_Iso193Rep2_buc,AIC_Iso193Rep2_gom,AIC_Iso193Rep2_bar)
rownames(AIC)<-c("Iso193Rep1_buccoef","Iso193Rep1_gomcoef","Iso193Rep1_barcoef",
                 "Iso193Rep2_buccoef","Iso193Rep2_gomcoef","Iso193Rep2_barcoef")
colnames(AIC)<-c("AIC")
BIC<-rbind(BIC_Iso193Rep1_buc,BIC_Iso193Rep1_gom,BIC_Iso193Rep1_bar,
           BIC_Iso193Rep2_buc,BIC_Iso193Rep2_gom,BIC_Iso193Rep2_bar)
rownames(BIC)<-c("Iso193Rep1_buccoef","Iso193Rep1_gomcoef","Iso193Rep1_barcoef",
                 "Iso193Rep2_buccoef","Iso193Rep2_gomcoef","Iso193Rep2_barcoef")
colnames(BIC)<-c("BIC")
Result_193<-cbind(Iso193parameters,AIC,BIC)

# Combine RMSE for each isolate
RMSE_193<-rbind(Iso193rmse.buc,Iso193rmse.gom,Iso193rmse.bar)
colnames(RMSE_193)<-c("RMSE")

# Generate output
Result_10dC<-rbind(Result_193,Result_194,Result_402,Result_407,Result_413,Result_433,
              Result_457,Result_474,Result_495,Result_518,Result_536,Result_564,
              Result_570,Result_638,Result_649)
as.data.frame(Result_10dC)
write.csv(Result_10dC,"Result_10dC.csv")

RMSE_10dC<-rbind(RMSE_193,RMSE_194,RMSE_402,RMSE_407,RMSE_413,RMSE_433,RMSE_457,
            RMSE_474,RMSE_495,RMSE_518,RMSE_536,RMSE_564,RMSE_570,RMSE_638,RMSE_649)
as.data.frame(RMSE_10dC)
write.csv(RMSE_10dC,"RMSE_10dC.csv")

# Compute combined RMSE 
resid.22.buc = sum(Iso193resid.buc^2,Iso194resid.buc^2,Iso402resid.buc^2,Iso407resid.buc^2,Iso413resid.buc^2,
                   Iso433resid.buc^2,Iso457resid.buc^2,Iso474resid.buc^2,Iso495resid.buc^2,Iso518resid.buc^2,
                   Iso536resid.buc^2,Iso564resid.buc^2,Iso570resid.buc^2,Iso638resid.buc^2,Iso649resid.buc^2)
length.22.buc = length(Iso193resid.buc)+length(Iso194resid.buc)+length(Iso402resid.buc)+length(Iso407resid.buc)+
                length(Iso413resid.buc)+length(Iso433resid.buc)+length(Iso457resid.buc)+length(Iso474resid.buc)+
                length(Iso495resid.buc)+length(Iso518resid.buc)+length(Iso536resid.buc)+length(Iso564resid.buc)+
                length(Iso570resid.buc)+length(Iso638resid.buc)+length(Iso649resid.buc)

resid.gom_tot = sum(resid.10.gom,resid.22.gom)
length.gom_tot = length.10.gom + length.22.gom
rmse.gom_tot = sqrt(resid.gom_tot/length.gom_tot)
RMSE_combined = c(rmse.bar_tot,rmse.buc_tot,rmse.gom_tot)

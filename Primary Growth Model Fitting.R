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
data_10dC$t<-as.numeric(data_10dC$t)
data_10dC$LOG10N<-as.numeric(data_10dC$LOG10N)


# Load primary growth models
source("UtilityFunctions.R")


#Note: This is an example of Isolate457Rep2. The code applies to other samples by changing sample names

# Subset data 
Iso457Rep2 <- subset(data_10dC, Isolate == "457" & Rep =="rep2")

plot(Iso457Rep2$t,Iso457Rep2$LOG10N)
mod <- lm(Iso457Rep2$LOG10N ~ Iso457Rep2$t)
slope <- coef(mod) [2]
slope


# Fit primary growth models 
# Fit Buchanan 
Iso457Rep2.buc_LM<- nlsLM(LOG10N ~ buchanan_log10N(t,lag,mumax,LOG10N0,LOG10Nmax), data=Iso457Rep2,
                          start=list (
                          LOG10N0 = 3,
                          lag = 100, 
                          mumax = 0.01, 
                          LOG10Nmax = 6), 
                          lower = c(0,0,0,0))

# Fit Gompertz 
Iso457Rep2.gom_LM <- nlsLM(LOG10N ~ gompertz_log10N(t,lag,mumax,LOG10N0,LOG10Nmax), data=Iso457Rep2,
                           start=list (
                           LOG10N0 = 3,
                           lag = 100,
                           mumax = 0.01,
                           LOG10Nmax = 6), 
                           lower = c(0,0,0,0))

# Fit Baranyi 
Iso457Rep2.bar_LM <- nlsLM(LOG10N ~ baranyi_log10N(t,lag,mumax,LOG10N0,LOG10Nmax), data=Iso457Rep2,
                           start=list (
                           LOG10N0 = 3,
                           lag = 100, 
                           mumax = 0.01, 
                           LOG10Nmax = 6), 
                           lower = c(0,0,0,0))

# Collect growth parameters estimated from 3 primary growth models
Iso457Rep2_buccoef<-coef(Iso457Rep2.buc_LM)
Iso457Rep2_gomcoef<-coef(Iso457Rep2.gom_LM)
Iso457Rep2_barcoef<-coef(Iso457Rep2.bar_LM)

# Compute AIC, BIC, RMSE for 3 primary growth models
AIC_Iso457Rep2_buc=AIC(Iso457Rep2.buc_LM)
BIC_Iso457Rep2_buc=BIC(Iso457Rep2.buc_LM)
Iso457resid.buc = c(Iso457Rep1.buc_LM$m$resid(), Iso457Rep2.buc_LM$m$resid())
Iso457rmse.buc = sqrt(sum(Iso457resid.buc^2)/length(Iso457resid.buc))

AIC_Iso457Rep2_gom=AIC(Iso457Rep2.gom_LM)
BIC_Iso457Rep2_gom=BIC(Iso457Rep2.gom_LM)
Iso457resid.gom = c(Iso457Rep1.gom_LM$m$resid(), Iso457Rep2.gom_LM$m$resid())
Iso457rmse.gom = sqrt(sum(Iso457resid.gom^2)/length(Iso457resid.gom))


AIC_Iso457Rep2_bar=AIC(Iso457Rep2.bar_LM)
BIC_Iso457Rep2_bar=BIC(Iso457Rep2.bar_LM)
Iso457resid.bar = c(Iso457Rep1.bar_LM$m$resid(), Iso457Rep2.bar_LM$m$resid())
Iso457rmse.bar = sqrt(sum(Iso457resid.bar^2)/length(Iso457resid.bar))

#Note: Isolate193Rep1&2 has an estimated lag time of 0 after fitting into Buchanan, Gompertz, Baranyi models
#Isolate193Rep1&2 are then fitted into Buchanan_without-lag and Baranyi_without_lag
#Isolate193 has to be repeated 

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
Iso495parameters<- rbind(Iso495Rep1_buccoef,Iso495Rep1_gomcoef,Iso495Rep1_barcoef,
                         Iso495Rep2_buccoef,Iso495Rep2_gomcoef,Iso495Rep2_barcoef)
AIC<-rbind(AIC_Iso495Rep1_buc,AIC_Iso495Rep1_gom,AIC_Iso495Rep1_bar,
           AIC_Iso495Rep2_buc,AIC_Iso495Rep2_gom,AIC_Iso495Rep2_bar)
rownames(AIC)<-c("Iso495Rep1_buccoef","Iso495Rep1_gomcoef","Iso495Rep1_barcoef",
                 "Iso495Rep2_buccoef","Iso495Rep2_gomcoef","Iso495Rep2_barcoef")
colnames(AIC)<-c("AIC")
BIC<-rbind(BIC_Iso495Rep1_buc,BIC_Iso495Rep1_gom,BIC_Iso495Rep1_bar,
           BIC_Iso495Rep2_buc,BIC_Iso495Rep2_gom,BIC_Iso495Rep2_bar)
rownames(BIC)<-c("Iso495Rep1_buccoef","Iso495Rep1_gomcoef","Iso495Rep1_barcoef",
                 "Iso495Rep2_buccoef","Iso495Rep2_gomcoef","Iso495Rep2_barcoef")
colnames(BIC)<-c("BIC")
Result_495<-cbind(Iso495parameters,AIC,BIC)

# Combine RMSE for each isolate
RMSE_495<-rbind(Iso495rmse.buc,Iso495rmse.gom,Iso495rmse.bar)
colnames(RMSE_495)<-c("RMSE")

#Generate output
Result_10dC<-rbind(Result_193,Result_194,Result_402,Result_407,Result_413,Result_433,
              Result_457,Result_474,Result_495,Result_518,Result_536,Result_564,
              Result_570,Result_638,Result_649)
as.data.frame(Result_10dC)
write.csv(Result_10dC,"Result_10dC.csv")

RMSE_10dC<-rbind(RMSE_193,RMSE_194,RMSE_402,RMSE_407,RMSE_413,RMSE_433,RMSE_457,
            RMSE_474,RMSE_495,RMSE_518,RMSE_536,RMSE_564,RMSE_570,RMSE_638,RMSE_649)
as.data.frame(RMSE_10dC)
write.csv(RMSE_10dC,"RMSE_10dC.csv")




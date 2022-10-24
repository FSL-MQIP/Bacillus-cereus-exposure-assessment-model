library(tidyverse)
library(stats)
install.packages('qpcR')
library('qpcR')

# Compute AIC and Akaike Weight

AIC_Iso433Rep2_buc_nolag=AIC(Iso433Rep2.buc_nolag_LM)
AIC_Iso433Rep2_bar_nolag=AIC(Iso433Rep2.bar_nolag_LM)
AIC_Iso433Rep2_nolag<-c(AIC_Iso433Rep2_buc_nolag,AIC_Iso433Rep2_bar_nolag)
Akaike_Weight_433Rep2_nolag<-akaike.weights(AIC_Iso433Rep2_nolag)$weights

Akaike_Weight<-cbind(Akaike_Weight_193Rep1,Akaike_Weight_193Rep2,Akaike_Weight_194Rep1,Akaike_Weight_194Rep2,
                     Akaike_Weight_125Rep1,Akaike_Weight_125Rep2,Akaike_Weight_135Rep1,Akaike_Weight_135Rep2,
                     Akaike_Weight_402Rep1,Akaike_Weight_402Rep2,Akaike_Weight_407Rep1,Akaike_Weight_407Rep2,
                     Akaike_Weight_413Rep1,Akaike_Weight_413Rep2,Akaike_Weight_433Rep1,Akaike_Weight_433Rep2,
                     Akaike_Weight_457Rep1,Akaike_Weight_457Rep2,Akaike_Weight_474Rep1,Akaike_Weight_474Rep2,
                     Akaike_Weight_495Rep1,Akaike_Weight_495Rep2,Akaike_Weight_518Rep1,Akaike_Weight_518Rep2,
                     Akaike_Weight_536Rep1,Akaike_Weight_536Rep2,Akaike_Weight_564Rep1,Akaike_Weight_564Rep2,
                     Akaike_Weight_570Rep1,Akaike_Weight_570Rep2,Akaike_Weight_638Rep1,Akaike_Weight_638Rep2,
                     Akaike_Weight_649Rep1,Akaike_Weight_649Rep2)

Akaike_Weight_433_nolag<-cbind(Akaike_Weight_433Rep1_nolag,Akaike_Weight_433Rep2_nolag)

write.csv(Akaike_Weight,"Akaike_Weight.csv")
write.csv(Akaike_Weight_433_nolag,"Akaike_Weight_433_nolag.csv")

# Import Dataset 

Akaike_Weight <-read.csv("Akaike_Weight.csv")
Akaike_Weight_433_nolag<-read.csv("Akaike_Weight_433_nolag.csv")
Result1<-read.csv("Result1.csv")
Result_433_nolag<-read.csv("Result_433_nolag.csv")

# Reshape Data 

AW_649_Rep2=Akaike_Weight[,35]
AW_649_Rep2$isolate=649
AW_649_Rep2$rep=2
AW_649_Rep2<-as.data.frame(AW_649_Rep2)
colnames(AW_649_Rep2)<-c('buc','gom','bar','Isolate','Rep')
AW<-rbind(AW_193_Rep1,AW_193_Rep2,AW_194_Rep1,AW_194_Rep2,
          AW_125_Rep1,AW_125_Rep2,AW_135_Rep1,AW_135_Rep2,
          AW_402_Rep1,AW_402_Rep2,AW_407_Rep1,AW_407_Rep2,
          AW_413_Rep1,AW_413_Rep2,AW_474_Rep1,AW_474_Rep2,
          AW_495_Rep1,AW_495_Rep2,AW_518_Rep1,AW_518_Rep2,
          AW_536_Rep1,AW_536_Rep2,AW_564_Rep1,AW_564_Rep2,
          AW_570_Rep1,AW_570_Rep2,AW_638_Rep1,AW_638_Rep2,
          AW_649_Rep1,AW_649_Rep2)
View(AW)

AW_433_Rep2=Akaike_Weight_433_nolag[,3]
AW_433_Rep2$isloate=433
AW_433_Rep2$rep=2
AW_433_Rep2<-as.data.frame(AW_433_Rep2)
colnames(AW_433_Rep2)<-c('buc','bar','Isolate','Rep')
AW_433<-rbind(AW_433_Rep1,AW_433_Rep2)

Result1<-Result1[,3:8]
colnames(Result1)<-c('lag','mumax','LOG10Nmax','Isolate','Rep','Model')
View(Result1)

Result_433_nolag<-Result_433_nolag[,3:8]
colnames(Result_433_nolag)<-c('lag','mumax','LOG10Nmax','Isolate','Rep','Model')

# Calculate weighted growth parameters 

Iso649Rep2_buc<- subset(Result1, Isolate == "649" & Rep =="2" & Model =="buc")
Iso649Rep2_gom<- subset(Result1, Isolate == "649" & Rep =="2" & Model =="gom")
Iso649Rep2_bar<- subset(Result1, Isolate == "649" & Rep =="2" & Model =="bar")
AW649Rep2<- subset(AW, Isolate == "649" & Rep =="2")

Iso649_lag<-Iso649Rep2_buc$lag*AW649Rep2$buc+Iso649Rep2_gom$lag*AW649Rep2$gom+Iso649Rep2_bar$lag*AW649Rep2$bar
Iso649Rep2_mumax<-Iso649Rep2_buc$mumax*AW649Rep2$buc+Iso649Rep2_gom$mumax*AW649Rep2$gom+Iso649Rep2_bar$mumax*AW649Rep2$bar
Iso649Rep2_LOG10Nmax<-Iso649Rep2_buc$LOG10Nmax*AW649Rep2$buc+Iso649Rep2_gom$LOG10Nmax*AW649Rep2$gom+Iso649Rep2_bar$LOG10Nmax*AW649Rep2$bar

Iso649Rep2gro<-cbind(Iso649Rep2_lag,Iso649Rep2_mumax,Iso649Rep2_LOG10Nmax)

Iso433Rep2_buc<- subset(Result_433_nolag, Isolate == "433" & Rep =="2" & Model =="buc")
Iso433Rep2_bar<- subset(Result_433_nolag, Isolate == "433" & Rep =="2" & Model =="bar")
AW433Rep2<- subset(AW_433, Isolate == "433" & Rep =="2")

Iso433Rep2_lag<-Iso433Rep2_buc$lag*AW433Rep2$buc+Iso433Rep2_bar$lag*AW433Rep2$bar
Iso433Rep2_mumax<-Iso433Rep2_buc$mumax*AW433Rep2$buc+Iso433Rep2_bar$mumax*AW433Rep2$bar
Iso433Rep2_LOG10Nmax<-Iso433Rep2_buc$LOG10Nmax*AW433Rep2$buc+Iso433Rep2_bar$LOG10Nmax*AW433Rep2$bar

Iso433Rep2gro<-cbind(Iso433Rep2_lag,Iso433Rep2_mumax,Iso433Rep2_LOG10Nmax)


Weighted_Growth<-cbind(Iso193Rep1gro,Iso194Rep1gro,Iso125Rep1gro,Iso135Rep1gro,Iso402Rep1gro,Iso407Rep1gro,
                       Iso413Rep1gro,Iso433Rep1gro,Iso474Rep1gro,Iso495Rep1gro,Iso518Rep1gro,Iso536Rep1gro,
                       Iso564Rep1gro,Iso570Rep1gro,Iso638Rep1gro,Iso649Rep1gro,Iso193Rep2gro,Iso194Rep2gro,
                       Iso125Rep2gro,Iso135Rep2gro,Iso402Rep2gro,Iso407Rep2gro,Iso413Rep2gro,Iso433Rep2gro,
                       Iso474Rep2gro,Iso495Rep2gro,Iso518Rep2gro,Iso536Rep2gro,Iso564Rep2gro,Iso570Rep2gro,
                       Iso638Rep2gro,Iso649Rep2gro)

Weighted_Growth<-as.data.frame(Weighted_Growth)
write.csv(Weighted_Growth,"Weighted_Growth")
View(Weighted_Growth)

Weighted_Growth_649Rep2=Weighted_Growth[,94:96]
Weighted_Growth_649Rep2$isloate=649
Weighted_Growth_649Rep2$rep=2
colnames(Weighted_Growth_649Rep2)<-c('lag','mumax','LOG10Nmax','Isolate','Rep')
rbind(Weighted_Growth_193Rep1,Weighted_Growth_193Rep2)

Weighted_Growth<-rbind(Weighted_Growth_193Rep1,Weighted_Growth_194Rep1,Weighted_Growth_125Rep1,Weighted_Growth_135Rep1,
                       Weighted_Growth_402Rep1,Weighted_Growth_407Rep1,Weighted_Growth_413Rep1,Weighted_Growth_433Rep1,
                       Weighted_Growth_474Rep1,Weighted_Growth_495Rep1,Weighted_Growth_518Rep1,Weighted_Growth_536Rep1,
                       Weighted_Growth_564Rep1,Weighted_Growth_570Rep1,Weighted_Growth_638Rep1,Weighted_Growth_649Rep1, 
                       Weighted_Growth_193Rep2,Weighted_Growth_194Rep2,Weighted_Growth_125Rep2,Weighted_Growth_135Rep2,
                       Weighted_Growth_402Rep2,Weighted_Growth_407Rep2,Weighted_Growth_413Rep2,Weighted_Growth_433Rep2,
                       Weighted_Growth_474Rep2,Weighted_Growth_495Rep2,Weighted_Growth_518Rep2,Weighted_Growth_536Rep2,
                       Weighted_Growth_564Rep2,Weighted_Growth_570Rep2,Weighted_Growth_638Rep2,Weighted_Growth_649Rep2)

write.csv(Weighted_Growth,"Weighted_Growth.csv")


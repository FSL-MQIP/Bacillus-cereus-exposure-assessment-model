library(tidyverse)
library('qpcR')

# Compute Akaike Weight
# Note: This is an example of Isolate536Rep3. The code applies to other samples by changing sample names
AIC_Iso536Rep3<-c(AIC_Iso536Rep3_buc,AIC_Iso536Rep3_gom,AIC_Iso536Rep3_bar)
Akaike_Weight_536Rep3<-akaike.weights(AIC_Iso536Rep3)$weights
Akaike_Weight_457Rep1<-as.data.frame(Akaike_Weight_457Rep1)
rownames(Akaike_Weight_457Rep1)=c("buc","gom","bar")

# Note: Result.csv was manually revised to add columns "Isolate","Rep","Model" in Excel and saved as "Result1.csv"
# Import dataset
Result1<-read.csv("Result1.csv")
Result1<-Result1[,3:8]
colnames(Result1)<-c('lag','mumax','LOG10Nmax','Isolate','Rep','Model')
View(Result1)

# Calculate weighted growth parameters 
Iso536Rep3_buc<- subset(Result1, Isolate == "536" & Rep =="3" & Model =="buc")
Iso536Rep3_gom<- subset(Result1, Isolate == "536" & Rep =="3" & Model =="gom")
Iso536Rep3_bar<- subset(Result1, Isolate == "536" & Rep =="3" & Model =="bar")

Iso536Rep3_lag<-Iso536Rep3_buc$lag*Akaike_Weight_536Rep3[1,1]+Iso536Rep3_gom$lag*Akaike_Weight_536Rep3[2,1]+Iso536Rep3_bar$lag*Akaike_Weight_536Rep3[3,1]
Iso536Rep3_mumax<-Iso536Rep3_buc$mumax*Akaike_Weight_536Rep3[1,1]+Iso536Rep3_gom$mumax*Akaike_Weight_536Rep3[2,1]+Iso536Rep3_bar$mumax*Akaike_Weight_536Rep3[3,1]
Iso536Rep3_LOG10Nmax<-Iso536Rep3_buc$LOG10Nmax*Akaike_Weight_536Rep3[1,1]+Iso536Rep3_gom$LOG10Nmax*Akaike_Weight_536Rep3[2,1]+Iso536Rep3_bar$LOG10Nmax*Akaike_Weight_536Rep3[3,1]

Iso536Rep3gro<-cbind(Iso536Rep3_lag,Iso536Rep3_mumax,Iso536Rep3_LOG10Nmax)
colnames(Iso649Rep2gro)=c("lag","mumax","LOG10Nmax")
rownames(Iso649Rep2gro)=c("Iso649Rep2")


Weighted_Growth<-rbind(Iso193Rep1gro,Iso194Rep1gro,Iso125Rep1gro,Iso135Rep1gro,Iso402Rep1gro,Iso407Rep1gro,
                       Iso413Rep1gro,Iso433Rep1gro,Iso457Rep3gro,Iso474Rep1gro,Iso495Rep1gro,Iso518Rep1gro,
                       Iso536Rep2gro,Iso564Rep1gro,Iso570Rep1gro,Iso638Rep1gro,Iso649Rep1gro,Iso193Rep2gro,
                       Iso194Rep2gro,Iso125Rep2gro,Iso135Rep2gro,Iso402Rep2gro,Iso407Rep2gro,Iso413Rep2gro,
                       Iso433Rep2gro,Iso457Rep4gro,Iso474Rep2gro,Iso495Rep2gro,Iso518Rep2gro,Iso536Rep3gro,
                       Iso564Rep2gro,Iso570Rep2gro,Iso638Rep2gro,Iso649Rep2gro)

Weighted_Growth<-as.data.frame(Weighted_Growth)
write.csv(Weighted_Growth,"Weighted_Growth.csv")




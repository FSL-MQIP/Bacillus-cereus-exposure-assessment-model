library(tidyverse)

data <-read.csv("Bacillus cereus growth curves 22 C new.csv")

#Note: This is a example of Isolate536. The code applies to other samples by changing sample names and column numbers

isolate536=data[3:length(data[,1])-1,49:51]
names(isolate536)=c("time","rep1","rep2")
isolate536=isolate536[1:16,]
isolate536=gather(isolate536,key = "rep", value = "count",2:3)
isolate536$isolate=536

data1<-rbind(isolate125,isolate135,isolate193,isolate194,isolate402,
             isolate407,isolate413,isolate433,isolate457,isolate474,
             isolate495,isolate518,isolate536,isolate564,isolate570,
             isolate638,isolate649)
data1<-as.data.frame(data1)
data1$time<-as.numeric(data1$time)
data1$count<-as.numeric(data1$count)
data1 <- na.omit(data1)
colnames(data1)<-c('t','Rep','LOG10N','Isolate')
write.csv(data1,"data1.csv")

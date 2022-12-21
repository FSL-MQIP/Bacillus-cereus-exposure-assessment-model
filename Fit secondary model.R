data <- read.csv("mumax.csv")
colnames(data) <- c("Isolate","T","mumax","sqrt_mumax")

Iso413<-subset(data,Isolate == "413")
fit <- lm(Iso413$sqrt_mumax[2:4] ~ Iso413$T[2:4])
plot(Iso413$T,Iso413$sqrt_mumax,
     ylim = c(0,1.5),
     xlim = c(0,24))
abline(fit)
summary(fit)

sqrt_mumax = 0.020914*T + 0.145254
  
b_413 <- coef(fit) [2]
Tmin_413 <- -0.145254/0.020914

data <- read.csv("mumax.csv")
colnames(data) <- c("Isolate","T","mumax","sqrt_mumax")

Iso193<-subset(data,Isolate == "193")
fit <- lm(Iso193$sqrt_mumax ~ Iso193$T)
plot(Iso193$T,Iso193$sqrt_mumax,
     ylim = c(0,1.5),
     xlim = c(0,24))
abline(fit)
summary(fit)

sqrt_mumax = 0.04495 *T - 0.20936
  
b_193 <- coef(fit) [2]
Tmin_193 <- 0.20936/0.04495

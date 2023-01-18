data <- read.csv("mumax_R.csv")
colnames(data) <- c("Isolate","T","mumax","sqrt_mumax")

Iso649<-subset(data,Isolate == "649")
fit <- lm(Iso649$sqrt_mumax ~ Iso649$T)
plot(Iso457$T,Iso457$sqrt_mumax,
     ylim = c(0,4),
     xlim = c(0,24))
abline(fit)
summary(fit)

sqrt_mumax = 0.067519 *T + 0.468946
  
b_649 <- coef(fit) [2]
Tmin_407 <- -0.468946/0.067519

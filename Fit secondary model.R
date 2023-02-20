data <- read.csv("mumax_R.csv")
colnames(data) <- c("Isolate","T","mumax","sqrt_mumax")

Iso649<-subset(data,Isolate == "649")
fit <- lm(Iso649$sqrt_mumax ~ Iso649$T)
plot(Iso649$T,Iso649$sqrt_mumax,
     ylim = c(0,5),
     xlim = c(0,24),
     main = "Isolate 649",
     xlab = "Temperature, °C",
     ylab = "Square root of mumax,sqrt(log CFU/mL per day)")
abline(fit)
summary(fit)

sqrt_mumax = 0.13263 *T - 0.8159
  
b_649 <- coef(fit) [2]
Tmin_649 <- 0.8159/0.13263

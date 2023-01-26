data <- read.csv("mumax_R.csv")
colnames(data) <- c("Isolate","T","mumax","sqrt_mumax")

Iso194<-subset(data,Isolate == "194")
fit <- lm(Iso194$sqrt_mumax ~ Iso194$T)
plot(Iso194$T,Iso194$sqrt_mumax,
     ylim = c(0,5),
     xlim = c(0,24),
     main = "Isolate 194 Secondary Model",
     xlab = "Temperature, °C",
     ylab = "Square root of mumax, ln CFU/mL h-1)1/2")
abline(fit)
summary(fit)

sqrt_mumax = 0.12459 *T - 0.65114
  
b_649 <- coef(fit) [2]
Tmin_649 <- 0.65114/0.12459

# Note: Weighted_Average vs Best_Model.csv was manually created by combining the results from Result1.csv, RMSE.csv and Weighted_Growth.csv in Excel
# The best model was selected by the smallest RMSE

# Import dataset
dat_wt<-read.csv("Weighted_Average vs Best_Model.csv")

# Association test: lag_wt vs lag_best
plot(dat_wt$lag_wt,dat_wt$lag_best,
     xlim = c(0,25),
     ylim = c(0,25),
     xlab="lag_best",
     ylab="lag_wt",
     main="lag")
reg<-lm(dat_wt$lag_wt ~ dat_wt$lag_best)
abline(reg,col="red")
m1=lm(dat_wt$lag_wt ~ dat_wt$lag_best)
slope_lag=coef(m1)[2]
cor_lag=cor(dat_wt$lag_wt,dat_wt$lag_best)

# Association test: mumax_wt vs mumax_best
plot(dat_wt$mumax_wt,dat_wt$mumax_best,
     xlim = c(0,2),
     ylim = c(0,2),
     xlab="mumax_best",
     ylab="mumax_wt",
     main="mumax")
reg<-lm(dat_wt$mumax_wt ~ dat_wt$mumax_best)
abline(reg,col="red")
m2=lm(dat_wt$mumax_wt ~ dat_wt$mumax_best)
slope_mumax=coef(m2)[2]
cor_mumax=cor(dat_wt$lag_wt,dat_wt$lag_best)

# Association test: LOG10Nmax_wt vs LOG10Nmax_best
plot(dat_wt$LOG10Nmax_wt,dat_wt$LOG10Nmax_best,
     xlim = c(0,10),
     ylim = c(0,10),
     xlab="LOG10Nmax_best",
     ylab="LOG10Nmax_wt",
     main="LOG10Nmax")
reg<-lm(dat_wt$LOG10Nmax_wt ~ dat_wt$LOG10Nmax_best)
abline(reg,col="red")
m3=lm(dat_wt$LOG10Nmax_wt ~ dat_wt$LOG10Nmax_best)
slope_LOG10Nmax=coef(m3)[2]
cor_LOG10Nmax=cor(dat_wt$LOG10Nmax_wt,dat_wt$LOG10Nmax_best)

correlation_coeff=cbind(c(cor_lag,cor_mumax,cor_LOG10Nmax))
slope=cbind(c(slope_lag,slope_mumax,slope_LOG10Nmax))
rownames(correlation_coeff)<-c("lag","mumax","LOG10Nmax")
colnames(correlation_coeff)<-c("cor")
rownames(slope)<-c("lag","mumax","LOG10Nmax")
colnames(slope)<-c("slope")
asso_test=cbind(correlation_coeff,slope)
write.csv(asso_test,"correlation_coeff.csv")



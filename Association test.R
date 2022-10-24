# Note: Weighted_Average vs Best_Model.csv was manually created by combining the results from Result1.csv, RMSE.csv and Weighted_Growth.csv in Excel
# The best model was selected by the smallest RMSE

# Import dataset
dat_wt<-read.csv("Weighted_Average vs Best_Model.csv")

# Association test: lag_wt vs lag_best
plot(dat_wt$lag_wt,dat_wt$lag_best,
     xlab="lag_best",
     ylab="lag_wt",
     main="lag")
reg<-lm(dat_wt$lag_wt ~ dat_wt$lag_best)
abline(reg,col="red")
cor_lag=cor(dat_wt$lag_wt,dat_wt$lag_best)

# Association test: mumax_wt vs mumax_best
plot(dat_wt$mumax_wt,dat_wt$mumax_best,
     xlab="mumax_best",
     ylab="mumax_wt",
     main="mumax")
reg<-lm(dat_wt$mumax_wt ~ dat_wt$mumax_best)
abline(reg,col="red")
cor_mumax=cor(dat_wt$lag_wt,dat_wt$lag_best)

# Association test: LOG10Nmax_wt vs LOG10Nmax_best
plot(dat_wt$LOG10Nmax_wt,dat_wt$LOG10Nmax_best,
     xlab="LOG10Nmax_best",
     ylab="LOG10Nmax_wt",
     main="LOG10Nmax")
reg<-lm(dat_wt$LOG10Nmax_wt ~ dat_wt$LOG10Nmax_best)
abline(reg,col="red")
cor_LOG10Nmax=cor(dat_wt$LOG10Nmax_wt,dat_wt$LOG10Nmax_best)

correlation_coeff=cbind(c(cor_lag,cor_mumax,cor_LOG10Nmax))
rownames(correlation_coeff)<-c("lag","mumax","LOG10Nmax")
write.csv(correlation_coeff,"correlation_coeff.csv")



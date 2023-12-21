setwd("C:/Users/sujun/Documents/GitHub/Bacillus-cereus-exposure-assessment-model/Exposure assessment model/plot")

library(dplyr)
library(dunn.test)

dat <- read.csv("GP.csv")
dat <- dat[dat$Group != "III", ]

dat_10 <- subset(dat, Temp == "10")
dat_22 <- subset(dat, Temp == "22")

result_10 <- dat_10 %>%
  group_by(Group) %>%
  summarise(
    lag_min = min(lag),
    lag_q25 = quantile(lag, 0.25),
    lag_median = median(lag),
    lag_q75 = quantile(lag, 0.75),
    lag_max = max(lag),
    mumax_min = min(mumax),
    mumax_q25 = quantile(mumax, 0.25),
    mumax_median = median(mumax),
    mumax_q75 = quantile(mumax, 0.75),
    mumax_max = max(mumax)
  )

result_22 <- dat_22 %>%
  group_by(Group) %>%
  summarise(
    lag_min = min(lag),
    lag_q25 = quantile(lag, 0.25),
    lag_median = median(lag),
    lag_q75 = quantile(lag, 0.75),
    lag_max = max(lag),
    mumax_min = min(mumax),
    mumax_q25 = quantile(mumax, 0.25),
    mumax_median = median(mumax),
    mumax_q75 = quantile(mumax, 0.75),
    mumax_max = max(mumax)
  )

kruskal.test(lag ~ Group, data = dat_10)
kruskal.test(mumax ~ Group, data = dat_10)

kruskal.test(lag ~ Group, data = dat_22)
dunn.test(dat_22$lag, dat_22$Group, method = "bonferroni")

kruskal.test(mumax ~ Group, data = dat_22)
dunn.test(dat_22$mumax, dat_22$Group, method = "bonferroni")


## Susie Data Test - hypoxia

lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel","data.table",
         "tidyverse","rstan","susieR", "glmnet","rethinking","monomvn"), require, character.only=T)

## Import data
df <- fread("../../data/GRDO_GEE_HA_NHD_2021_01_30.csv")
colnames(df)

## Create data needed for susie function: https://stephenslab.github.io/susieR/reference/susie.html
X <- subset(df, select = c(slope_calc:NHD_PctHbWet2011Ws, DO_mgL_mean))
X <- X %>% select_if(is.numeric)
X[mapply(is.infinite, X)] <- NA
X <- na.omit(X)

X <- sample_n(X, 1000)

## standardize
X_std <- as.data.frame(apply(X, 2, function(y) (y - mean(y))/sd(y)))

y <- X_std$DO_mgL_mean
X_std <- X_std[,-which(names(X_std) %in% c("DO_mgL_mean","NHD_REACHCODE","HYRIV_ID","NHD_COMID"))]

## Run susie
res <- susie(as.matrix(X_std), y, L=10)
summary(res)


plot(y, predict(res))


susie_get_posterior_mean(res)










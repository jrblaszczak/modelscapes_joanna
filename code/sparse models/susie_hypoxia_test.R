## Susie Data Test - hypoxia

lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel","data.table",
         "tidyverse","rstan","susieR", "glmnet","rethinking","monomvn"), require, character.only=T)

## Import data
df <- fread("../../data/GRDO_GEE_HA_NHD_2021_01_30.csv")
colnames(df)

## Create data needed for susie function: https://stephenslab.github.io/susieR/reference/susie.html
X <- subset(df, select = c(slope_calc:NHD_PctHbWet2011Ws, Hyp_pr_sub2))
X <- X %>% select_if(is.numeric)
X[mapply(is.infinite, X)] <- NA
X <- na.omit(X)

## standardize
X_std <- as.data.frame(apply(X, 2, function(y) (y - mean(y))/sd(y)))

y <- X_std$Hyp_pr_sub2
X_std$Hyp_pr_sub2 <- NULL

## Run susie
res <- susie(X, y, L=10)



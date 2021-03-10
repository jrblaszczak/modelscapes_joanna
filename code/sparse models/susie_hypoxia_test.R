## Susie Data Test - hypoxia

lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel","data.table",
         "tidyverse","rstan","susieR", "glmnet","rethinking","monomvn"), require, character.only=T)

## Import data
df <- fread("../../data/GRDO_GEE_HA_NHD_2021_01_30.csv")
colnames(df)

## Create data needed for susie function: https://stephenslab.github.io/susieR/reference/susie.html
X <- subset(df, select = c(slope_calc:NHD_PctHbWet2011Ws, DO_mgL_mean))
X <- X %>% select_if(is.numeric)
X <- na.omit(X)

y <- X$DO_mgL_mean
X$DO_mgL_mean <- NULL

## Run susie
res <- susie(X, y, L=10)



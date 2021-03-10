## Testing sparse modeling approaches on simulated and existing data

# load packages
#devtools::install_github("stephenslab/susieR")

lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm",
         "susieR", "glmnet","rethinking","monomvn"), require, character.only=T)

## Susie vignette: https://stephenslab.github.io/susieR/articles/index.html

## Simple example
set.seed(1)
n    <- 100 # rows
p    <- 100 # columns
beta <- rep(0,p) # coefficients
beta[c(1,2,50,92)] <- 1 # set a subset of coefficients to be particularly high
X   <- matrix(rnorm(n*p),nrow=n,ncol=p) # create a nxp matrix
y   <- X %*% beta + rnorm(n) ## matrix multiplication between beta and X to generate observed responses, a vector of length n
res <- susie(X,y,L=10) ## susie function: https://stephenslab.github.io/susieR/reference/susie.html
plot(coef(res),pch = 20)
plot(y,predict(res),pch = 20)


## Andrew Siefart
# simulate data -----------------------------------------------------------

# simulate correlated predictor variables from multivariate normal distribution
n <- 1000
p <- 800
mu <- rnorm(p, 0, 1)                 # means
Sigma <- rethinking::rlkjcorr(1, p)  # covariance matrix
X <- MASS::mvrnorm(n, mu = mu, Sigma = Sigma)  # draw predictors from multivariate normal distribution  

# regression coefficients; most close to zero but some large values
beta <- rgamma(p, 0.01, 0.1) * sample(c(-1, 1), p, replace = T)
plot(ecdf(beta)) #empirical cumulative distribution function (ecdf) is a step function with jumps at observation values,
#where is the number of tied observations at that value.

# simulate response variable
y <- X %*% beta + rnorm(n)

# fit models and plot actual vs. estimated coefficients -----------------------

# linear regression
linreg <- lm(y~X)
plot(abs(beta), abs(coef(linreg)[-1])); abline(0, 1)

# susie
res <- susie(X, y, L=10)
plot(abs(beta), abs(coef(res)[-1])); abline(0, 1)

# lasso
lasso <- glmnet(X, y, alpha = 1, nlambda = 3)
plot(abs(beta), abs(coef(lasso)[-1,1])); abline(0, 1)
plot(abs(beta), abs(coef(lasso)[-1,2])); abline(0, 1)
plot(abs(beta), abs(coef(lasso)[-1,3])); abline(0, 1)

# ridge
ridge <- glmnet(X, y, alpha = 0, nlambda = 3)
plot(abs(beta), abs(coef(ridge)[-1,1])); abline(0, 1)
plot(abs(beta), abs(coef(ridge)[-1,2])); abline(0, 1)
plot(abs(beta), abs(coef(ridge)[-1,3])); abline(0, 1)


#################
## Chhaya Werner
###################
## code to visualize how the results change with sample size n,
## number of parameters p, and number of causal effects L. 

set.seed(1)

# function to test how well susie performs with different inputs
# output is a data frame comparing susie results to a linear model
test_susie <- function(n, p, L){
  
  # choose the causal parameters randomly
  beta <- rep(0,p)
  causal <- sample(1:p, size = L)
  beta[causal] <- 1 
  X <- matrix(rnorm(n*p), nrow=n, ncol=p)
  
  y <- X %*% beta + rnorm(n)
  
  ## do sparse modeling with susie
  res <- susie(X, y, L=10) # keeping susie's expected L at 10
  # summary(res)
  
  ## linear model using the causal parameters
  lm.result <- lm(y ~ X[,causal])
  # summary(lm.result)
  
  ## combine into a data frame to compare
  df <- data.frame(sus.pred = predict(res), lm.pred = predict(lm.result),
                   sample = n, parameters = p, causal = L)
  return(df)
}

## running our function over various n and p combinations
n.test <- c(20, 50, 200, 1000)
p.test <- c(10, 20, 50, 200, 800)
L <- 5
susie.predict <- data.frame()
for(n in n.test){
  for(p in p.test){
    susie.predict.new <- test_susie(n, p, L)
    susie.predict <- rbind(susie.predict, susie.predict.new)
  }
}

## visualize the output
ggplot(susie.predict, aes(x = lm.pred, y = sus.pred)) +
  facet_grid(sample ~ parameters) + 
  geom_point() +
  theme_bw() + 
  xlab('lm prediction') + 
  ylab('susie prediction') +
  ggtitle("L = 5, cols = parameters, rows = samples")


## what if L is larger?
L <- 10
susie.predict <- data.frame()
for(n in n.test){
  for(p in p.test){
    susie.predict.new <- test_susie(n, p, L)
    susie.predict <- rbind(susie.predict, susie.predict.new)
  }
}

ggplot(susie.predict, aes(x = lm.pred, y = sus.pred)) +
  facet_grid(sample ~ parameters) + 
  geom_point() +
  theme_bw() + 
  xlab('lm prediction') + 
  ylab('susie prediction') +
  ggtitle("L = 10, cols = parameters, rows = samples")




#################
## Alex Buerkle
#################
blasso.result<-blasso(X,y, RJ=TRUE,M=12,T=25000)
plot(blasso.result$lpost) # plot the model likelihood
# plot posterior inclusion probabilities for coefficients
plot(summary(blasso.result, burnin=10000)$bn0) 
# model is not getting a good sparse set and ends up against max of 12 coeff 
plot(blasso.result, burnin=10000, which='m')
# model does a good job of estimating the effects of the causal covariates
plot(colMeans(blasso.result$beta[10001:25000,]))


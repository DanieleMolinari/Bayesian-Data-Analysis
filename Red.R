red <- read.csv("winequality-red.csv", header = TRUE, sep = ";")
red <- as.data.frame(red)
summary(red$fixed.acidity)

head(red)

hist(red$quality, xlab = "", main = "Red wine quality", ylab = "Quality", breaks = 5, xlim = c(0, 10), col = "#4575b4")

correlred <- cor(red)
corrplot::corrplot(correlred, method = "number", type = "upper")


#############################################################################
################Lasso with Metropolis - Hastings sampler#####################

Xred <- red[,-12]
yred <- red$quality
Xred <- matrix(unlist(Xred), ncol = dim(Xred)[2], byrow = FALSE)
Xred <- cbind(rep(1, nrow(Xred)), scale(Xred))
n <- length(yred)
k <- ncol(Xred)
npars <- k + (k-1) + 2
nchains <- 4

logprior <- function(pars) {
  b <- as.matrix(pars[1:k])
  tau <- pars[(k+1):(2*k-1)]
  lambda2 <- pars[2*k]
  sigma2 <- pars[npars]
  if (lambda2 <= 0)
    return(log(0))
  if (sigma2 <= 0)
    return(log(0))
  if (any(tau<0))
    return(log(0))
  dens <- (sum(dnorm(abs(b[2:k]),0,sqrt(sigma2*tau),log=TRUE))+
      dgamma(1/sigma2,1,5,log=TRUE) +
      dgamma(lambda2,0.02,0.1,log=TRUE)+
      sum(dexp(tau,lambda2/2,log=TRUE)))
  return(dens)
}

loglikelihood <- function(yred,Xred,par) {
  b <- as.matrix(par[1:k])
  sigma2 <- par[npars]
  if (sigma2 <= 0)
    return(log(0))
  dens <- sum(dnorm(yred,t(Xred%*%b),sigma2,log=TRUE))
  return(dens)
}

draws <- matrix(0,nrow=40000,ncol=npars)
draws[1,npars] <- 1/rgamma(1,1,15)
draws[1,npars-1] <- 1
draws[1,1] <- mean(yred)
draws[1,2:k] <- rep(0,k-1)
draws[1,(k+1):(2*k-1)] <- rep(1,k-1)
for (step in 2:40000) {
  proposed <- rnorm(npars,draws[step-1,],0.1)
  r <- loglikelihood(yred,Xred,proposed)+
    logprior(proposed)-
    loglikelihood(yred,Xred,draws[step-1,])-
    logprior(draws[step-1,])
  u <- runif(1)
  if (log(u) < r) {
    draws[step,] <- proposed
  } else {
    draws[step,] <- draws[step-1,]
  }
}
sample <- draws[20000:40000,]
betas <- sample[, 2:k]

boxplot(betas, outline = FALSE, horizontal = TRUE)
abline(v = 0, lty = 2)
#############################################################################
#############################################################################

library(rstanarm)
library(bayesplot)

#Model1
set.seed(10)
Red_stan_model <- stan_glm(quality ~ ., data = red, iter = 1000, warmup = 500)
summary(Red_stan_model)

dimnames(as.array(Red_stan_model))

color_scheme_set("blue")
mcmc_trace(as.array(Red_stan_model), pars = c("(Intercept)", "sigma"),
           n_warmup = 500)

red_Rhats <- rhat(Red_stan_model)
hat1 <- mcmc_rhat(red_Rhats)

#model 2
set.seed(10)
Red_stan_model2 <- stan_glm(quality ~ ., data = red, iter = 1000, warmup = 500, prior_intercept = normal(0, 5), prior = normal(0, 2))
summary(Red_stan_model2)

mcmc_trace(as.array(Red_stan_model2), pars = c("(Intercept)", "sigma"),
           n_warmup = 500)

red_Rhats2 <- rhat(Red_stan_model2)
hat2 <- mcmc_rhat(red_Rhats2)

#model 3
set.seed(10)
Red_stan_model3 <- stan_glm(quality ~ ., data = red, iter = 1000, warmup = 500, prior_intercept = normal(0, 15), prior = normal(0, 5))
summary(Red_stan_model3)

mcmc_trace(as.array(Red_stan_model3), pars = c("(Intercept)", "sigma"),
           n_warmup = 500)

red_Rhats3 <- rhat(Red_stan_model3)
hat3 <- mcmc_rhat(red_Rhats3)

#model 4
set.seed(10)
Red_stan_model4 <- stan_glm(quality ~ ., data = red, iter = 1000, warmup = 500, prior_intercept = normal(1, 10), prior = normal(1, 2.5))
summary(Red_stan_model4)


mcmc_trace(as.array(Red_stan_model4), pars = c("(Intercept)", "sigma"),
           n_warmup = 500)

red_Rhats4 <- rhat(Red_stan_model4)
hat4 <- mcmc_rhat(red_Rhats4)

#model 5
set.seed(10)
Red_stan_model5 <- stan_glm(quality ~ ., data = red, iter = 1000, warmup = 500, prior_intercept = normal(1, 5), prior = normal(1, 2))
summary(Red_stan_model5)

mcmc_trace(as.array(Red_stan_model5), pars = c("(Intercept)", "sigma"),
           n_warmup = 500)

red_Rhats5 <- rhat(Red_stan_model5)
hat5 <- mcmc_rhat(red_Rhats5)

#model 6
set.seed(10)
Red_stan_model6 <- stan_glm(quality ~ ., data = red, iter = 1000, warmup = 500, prior_intercept = normal(1, 15), prior = normal(1, 5))
summary(Red_stan_model6)
prior_summary(Red_stan_model6)
mcmc_trace(as.array(Red_stan_model6), pars = c("(Intercept)", "fixed.acidity", "volatile.acidity", "citric.acid",
                                               "residual.sugar", "chlorides", "free.sulfur.dioxide" ,
                                               "total.sulfur.dioxide", "density", "pH", "sulphates", "alcohol", "sigma"), n_warmup = 500)

red_Rhats6 <- rhat(Red_stan_model6)
hat6 <- mcmc_rhat(red_Rhats6)

library(cowplot)
plot_grid(hat1, hat2, hat3, hat4, hat5, hat6, ncol = 3, labels = c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6"), hjust = -1.5)


plot(Red_stan_model6, "acf", pars = c("(Intercept)", "fixed.acidity", "volatile.acidity", "citric.acid",
                                     "residual.sugar", "chlorides", "free.sulfur.dioxide" ,
                                     "total.sulfur.dioxide", "density", "pH", "sulphates", "alcohol", "sigma"))

red_cred_int <- posterior_interval(Red_stan_model6, prob = 0.95)
red_cred_int
mcmc_intervals(as.array(Red_stan_model6), pars = c("(Intercept)", "fixed.acidity", "volatile.acidity", "citric.acid",
                                                  "residual.sugar", "chlorides", "free.sulfur.dioxide" ,
                                                  "total.sulfur.dioxide", "density", "pH", "sulphates", "alcohol", "sigma"))
mcmc_intervals(as.array(Red_stan_model6), pars = c("volatile.acidity", "citric.acid",
                                                  "chlorides", "pH", "sulphates"))
mcmc_intervals(as.array(Red_stan_model6), pars = c("fixed.acidity", "residual.sugar", "free.sulfur.dioxide",
                                                  "total.sulfur.dioxide", "alcohol", "sigma"))
mcmc_intervals(as.array(Red_stan_model6), pars = c("free.sulfur.dioxide", "total.sulfur.dioxide"))

prior_summary(Red_stan_model6)
15*sd(red$quality)
5/(sd(red$fixed.acidity))*sd(red$quality)
5/(sd(red$volatile.acidity))*sd(red$quality)
5/(sd(red$citric.acid))*sd(red$quality)


#blasso uses Gibbs sampler
library(monomvn)

yred <- red$quality
Xred <- red[,-12]
Xred <- matrix(unlist(Xred), ncol = dim(Xred)[2], byrow = FALSE)
Xred <- scale(Xred)
redbla <- blasso(Xred, yred, T=15000, rd = c(1, 15), ab = c(1, 5),
                 thin = 100, RJ = FALSE, normalize = FALSE)

Red_plot <- boxplot(redbla$beta, outline = FALSE, horizontal = TRUE)
abline(v = 0, lty = 2) #1, 3, 4, 8 out

#selected model
red_fin <- red[, -c(1, 3, 4, 8)]
set.seed(10)
Red_stan_model_select <- stan_glm(quality ~ ., data = red_fin, iter = 1000, warmup = 500, prior_intercept = normal(1, 15), prior = normal(1, 5))
mcmc_trace(as.array(Red_stan_model_select), pars = c("(Intercept)", "volatile.acidity", 
                                                     "chlorides", "free.sulfur.dioxide", "total.sulfur.dioxide",
                                                     "pH", "sulphates", "alcohol", "sigma"), n_warmup = 500)

red_Rhats_select <- rhat(Red_stan_model_select)
hat_select <- mcmc_rhat(red_Rhats_select)
hat_select

plot(Red_stan_model_select, "acf", pars = c("(Intercept)", "volatile.acidity",
                                      "chlorides", "free.sulfur.dioxide", "total.sulfur.dioxide",
                                      "pH", "sulphates", "alcohol", "sigma"))


#predictions
predictions_firstmodel <- posterior_linpred(Red_stan_model6)

iteration_1_6 <- predictions_firstmodel[1,]
iteration_328_6 <- predictions_firstmodel[328,]
iteration_500_6 <- predictions_firstmodel[500,]
iteration_1025_6 <- predictions_firstmodel[1025,]


summary(red$quality)

summary(iteration_1_6)
summary(iteration_328_6)
summary(iteration_500_6)
summary(iteration_1025_6)

yrep_6 <- posterior_predict(Red_stan_model6)

par(mfrow = c(2, 3))

first_6 <- predictions_firstmodel[,1375]
summary(first_6)
red$quality[1375]
plot(density(yrep_6[,1375]), col = "red", main = "Quality Wine 3")
abline(v = red[1375, 12])

second_6 <- predictions_firstmodel[,928]
summary(second_6)
red[928,12]
plot(density(yrep_6[,928]), col = "red", main = "Quality Wine 4")
abline(v = red[928, 12])

third_6 <- predictions_firstmodel[,203]
summary(third_6)
red[203,12]
plot(density(yrep_6[,203]), col = "red", main = "Quality Wine 5")
abline(v = red[203, 12])

fourth_6 <- predictions_firstmodel[,238]
summary(fourth_6)
red[238,12]
plot(density(yrep_6[,238]), col = "red", main = "Quality Wine 6")
abline(v = red[238, 12])

fifth_6 <- predictions_firstmodel[,545]
summary(fifth_6)
red[545,12]
plot(density(yrep_6[,454]), col = "red", main = "Quality Wine 7")
abline(v = red[454, 12])

sixth_6 <- predictions_firstmodel[,1270]
summary(sixth_6)
red[1270,12]
plot(density(yrep_6[,1270]), col = "red", main = "Quality Wine 8")
abline(v = red[1270, 12])

par(mfrow = c(1,1))

pp_check(Red_stan_model6, "dens_overlay", adjust = 5)
pp_check(Red_stan_model6, "stat")
pp_check(Red_stan_model6, "stat_2d")


predictions_selectmodel <- posterior_linpred(Red_stan_model_select)

iteration_1_sel <- predictions_selectmodel[1,]
iteration_328_sel <- predictions_selectmodel[328,]
iteration_500_sel <- predictions_selectmodel[500,]
iteration_1025_sel <-predictions_selectmodel[1025,]
summary(iteration_1_sel)
summary(iteration_328_sel)
summary(iteration_500_sel)
summary(iteration_1025_sel)

yrep_sel <- posterior_predict(Red_stan_model_select)

par(mfrow = c(2, 3))

first_Sel <- predictions_selectmodel[,1375]
summary(first_Sel)
red$quality[1375]
plot(density(yrep_sel[,1325]), col = "red")
abline(v = red[1375, 12])

second_Sel <- predictions_selectmodel[,928]
summary(second_Sel)
red[928,12]
plot(density(yrep_sel[,928]), col = "red")
abline(v = red[928, 12])

third_Sel <- predictions_selectmodel[,203]
summary(third_Sel)
red[203,12]
plot(density(yrep_sel[,203]), col = "red")
abline(v = red[203, 12])

fourth_Sel <- predictions_selectmodel[,238]
summary(fourth_Sel)
red[238,12]
plot(density(yrep_sel[,238]), col = "red")
abline(v = red[238, 12])

fifth_Sel <- predictions_selectmodel[,545]
summary(fifth_Sel)
red[545,12]
plot(density(yrep_sel[,454]), col = "red")
abline(v = red[454, 12])

sixth_Sel <- predictions_selectmodel[,1270]
summary(sixth_Sel)
red[1270,12]
plot(density(yrep_sel[,1270]), col = "red")
abline(v = red[1270, 12])

par(mfrow = c(1, 1))

pp_check(Red_stan_model_select, "dens_overlay", adjust = 5)
pp_check(Red_stan_model_select, "stat")
pp_check(Red_stan_model_select, "stat_2d")

#Do the log first?
red$fixed.acidity[red$fixed.acidity==0] <- 0.1
red$volatile.acidity[red$volatile.acidity==0] <- 0.1
red$citric.acid[red$citric.acid==0] <- 0.1
red$residual.sugar[red$residual.sugar==0] <- 0.1
red$chlorides[red$chlorides==0] <- 0.1
red$free.sulfur.dioxide[red$free.sulfur.dioxide==0] <- 0.1
red$total.sulfur.dioxide[red$total.sulfur.dioxide==0] <- 0.1
red$density[red$density==0] <- 0.1
red$pH[red$pH==0] <- 0.1
red$sulphates[red$sulphates==0] <- 0.1
red$alcohol[red$alcohol==0] <- 0.1

log_Red_stan_model <- stan_glm(log(quality) ~ log(fixed.acidity) + log(red$volatile.acidity) + log(citric.acid) + 
                               log(residual.sugar) + log(chlorides) + log(free.sulfur.dioxide) + 
                               log(total.sulfur.dioxide) + log(density) + log(pH) + log(sulphates) +
                               log(alcohol), data = red, iter = 1000, warmup = 500, 
                               prior_intercept = normal(1, 15), prior = normal(1, 5), adapt_delta = 0.99)
summary(log_Red_stan_model)

mcmc_trace(as.array(log_Red_stan_model), pars = c("(Intercept)", "fixed.acidity", "volatile.acidity", "citric.acid",
                                                        "residual.sugar", "chlorides", "free.sulfur.dioxide" ,
                                                        "total.sulfur.dioxide", "density", "pH", "sulphates", "alcohol", "sigma"), n_warmup = 500)

red_Rhats_log_model <- rhat(log_Red_stan_model)
hat_log_model <- mcmc_rhat(red_Rhats_log_model)
hat_log_model

plot(log_Red_stan_model, "acf", pars = c("(Intercept)", "fixed.acidity", "volatile.acidity", "citric.acid",
                                         "residual.sugar", "chlorides", "free.sulfur.dioxide" ,
                                         "total.sulfur.dioxide", "density", "pH", "sulphates", "alcohol", "sigma"))

predictions_log_model <- posterior_linpred(log_Red_stan_model)

iteration_1_log <- predictions_log_model[1,]
iteration_328_log <- predictions_log_model[328,]
iteration_500_log <- predictions_log_model[500,]
iteration_1025_log <- predictions_log_model[1025,]

summary(log(red$quality))

summary(iteration_1_log)
summary(iteration_328_log)
summary(iteration_500_log)
summary(iteration_1025_log)

yrep_log <- posterior_predict(log_Red_stan_model)

par(mfrow = c(2, 3))

first_log <- predictions_log_model[,1375]
summary(first_log)
log(red$quality[1375])
plot(density(yrep_log[,1375]), col = "red",  main = "Quality Wine 3")
abline(v = log(red[1375, 12]))

second_log <- predictions_log_model[,928]
summary(second_log)
log(red[928,12])
plot(density(yrep_log[,928]), col = "red",  main = "Quality Wine 4")
abline(v = log(red[928, 12]))

third_log <- predictions_log_model[,203]
summary(third_log)
log(red[203,12])
plot(density(yrep_log[,203]), col = "red",  main = "Quality Wine 5")
abline(v = log(red[203, 12]))

fourth_log <- predictions_log_model[,238]
summary(fourth_log)
log(red[238,12])
plot(density(yrep_log[,238]), col = "red",  main = "Quality Wine 6")
abline(v = log(red[238, 12]))

fifth_log <- predictions_log_model[,545]
summary(fifth_log)
log(red[545,12])
plot(density(yrep_log[,454]), col = "red",  main = "Quality Wine 7")
abline(v = log(red[454, 12]))

sixth_log <- predictions_log_model[,1270]
summary(sixth_log)
log(red[1270,12])
plot(density(yrep_log[,1270]), col = "red",  main = "Quality Wine 8")
abline(v = log(red[1270, 12]))

par(mfrow = c(1, 1))

log_pp1 <- pp_check(log_Red_stan_model, "dens_overlay", adjust = 5)
log_pp2 <- pp_check(log_Red_stan_model, "stat")
log_pp3 <- pp_check(log_Red_stan_model, "stat_2d")

plot_grid(log_pp1, log_pp2, log_pp3, labels = c("a)", "b)", "c)"))

#Add pairwise variables
red$rspH <- red$residual.sugar*red$pH
red$catsd <- red$citric.acid*red$total.sulfur.dioxide
red$sa <- red$sulphates*red$alcohol
red$fac <- red$fixed.acidity/red$chlorides
red$dpH <- red$density/red$pH
red$fsda <- red$free.sulfur.dioxide/red$alcohol

#model transformed
red_stan_model_trans <- stan_glm(quality ~ ., data = red, iter = 1000, warmup = 500, prior_intercept = normal(1, 15), prior = normal(1, 5))
summary(red_stan_model_trans)

mcmc_trace(as.array(red_stan_model_trans), pars = c("(Intercept)", "fixed.acidity", "volatile.acidity", "citric.acid",
                                               "residual.sugar", "chlorides", "free.sulfur.dioxide" ,
                                               "total.sulfur.dioxide", "density", "pH", "sulphates", "alcohol", "rspH",
                                               "catsd", "sa", "fac", "dpH", "fsda", "sigma"), n_warmup = 500)

red_Rhats_model_trans <- rhat(red_stan_model_trans)

mcmc_rhat(red_Rhats_model_trans)

predictions_trans <- posterior_linpred(red_stan_model_trans)

yrep_trans <- posterior_predict(red_stan_model_trans)

par(mfrow = c(2, 3))

first_trans <- predictions_trans[,1375]
summary(first_trans)
red$quality[1375]
plot(density(yrep_trans[,1375]), col =  "red", main =  "Quality Wine 3")
abline(v = red[1375, 12])

second_trans <- predictions_trans[, 928]
summary(second_trans)
red[928,12]
plot(density(yrep_trans[,928]), col = "red", main =  "Quality Wine 4")
abline(v = red[928, 12])

third_trans <- predictions_trans[, 203]
summary(third_trans)
red[203, 12]
plot(density(yrep_trans[, 203]), col = "red", main =  "Quality Wine 5")
abline(v = red[203, 12])

fourth_trans <- predictions_trans[,238]
summary(fourth_trans)
red[238,12]
plot(density(yrep_trans[, 238]), col = "red", main = "Quality Wine 6")
abline(v = red[238, 12])

fifth_trans <- predictions_trans[, 545]
summary(fifth_trans)
red[545, 12]
plot(density(yrep_trans[, 545]), col = "red",  main = "Quality Wine 7")
abline(v = red[545, 12])

sixth_trans <- predictions_trans[, 1270]
summary(sixth_trans)
red[1270, 12]
plot(density(yrep_trans[, 1270]), col = "red",  main = "Quality Wine 8")
abline(v = red[1270, 12])

par(mfrow = c(1, 1))

pp_check(red_stan_model_trans, "dens_overlay", adjust = 5)
pp_check(red_stan_model_trans, "stat")
pp_check(red_stan_model_trans, "stat_2d")

#Reapplying lasso scratch...

#Reapplying Blasso function
library(monomvn)

yred <- red$quality
Xred <- red[,-12]
Xred <- matrix(unlist(Xred), ncol = dim(Xred)[2], byrow = FALSE)
Xred <- scale(Xred)
redbla <- blasso(Xred, yred, T=15000, rd = c(1, 15), ab = c(1, 5),
                 thin = 100, RJ = FALSE, normalize = FALSE)

Red_plot <- boxplot(redbla$beta, outline = FALSE, horizontal = TRUE)
abline(v = 0, lty = 2)

#taking away 1, 3, 8, 11

red_fin_2 <- red[, -c(1, 3, 8, 11)]

Red_stan_model_trans_sel <- stan_glm(quality ~ ., data = red_fin_2, iter = 1000, warmup = 500, prior_intercept = normal(1, 15), prior = normal(1, 5))

mcmc_trace(as.array(Red_stan_model_trans_sel), pars = c("(Intercept)", "volatile.acidity", 
                                                    "residual.sugar", "chlorides", "free.sulfur.dioxide" ,
                                                    "total.sulfur.dioxide", "pH", "sulphates", "rspH",
                                                    "catsd", "sa", "fac", "dpH", "fsda", "sigma"), n_warmup = 500)

red_Rhats_model_trans_sel <- rhat(Red_stan_model_trans_sel)

mcmc_rhat(red_Rhats_model_trans_sel)

predictions_trans_sel <- posterior_linpred(Red_stan_model_trans_sel)

yrep_trans_sel <- posterior_predict(Red_stan_model_trans_sel)

par(mfrow = c(2, 3))

first_trans_sel <- predictions_trans_sel[,1375]
summary(first_trans_sel)
red$quality[1375]
plot(density(yrep_trans_sel[,1325]), col = "red",  main =  "Quality Wine 3")
abline(v = red[1375, 12])

second_trans_sel <- predictions_trans_sel[,928]
summary(second_trans_sel)
red[928,12]
plot(density(yrep_trans_sel[,928]), col = "red",  main =  "Quality Wine 3")
abline(v = red[928, 12])

third_trans_sel <- predictions_trans_sel[,203]
summary(third_trans_sel)
red[203,12]
plot(density(yrep_trans_sel[,203]), col = "red",  main =  "Quality Wine 3")
abline(v = red[203, 12])

fourth_trans_sel <- predictions_trans_sel[,238]
summary(fourth_trans_sel)
red[238,12]
plot(density(yrep_trans_sel[,238]), col = "red",  main =  "Quality Wine 3")
abline(v = red[238, 12])

fifth_trans_sel <- predictions_trans_sel[,545]
summary(fifth_trans_sel)
red[545,12]
plot(density(yrep_trans_sel[,454]), col = "red",  main =  "Quality Wine 3")
abline(v = red[454, 12])

sixth_trans_sel <- predictions_trans_sel[,1270]
summary(sixth_trans_sel)
red[1270,12]
plot(density(yrep_trans_sel[,1270]), col = "red",  main =  "Quality Wine 3")
abline(v = red[1270, 12])

par(mfrow = c(1, 1))

pp_check(Red_stan_model_trans_sel, "dens_overlay", adjust = 5)
pp_check(Red_stan_model_trans_sel, "stat")
pp_check(Red_stan_model_trans_sel, "stat_2d")

#try log

red$fixed.acidity[red$fixed.acidity==0] <- 0.1
red$volatile.acidity[red$volatile.acidity==0] <- 0.1
red$citric.acid[red$citric.acid==0] <- 0.1
red$residual.sugar[red$residual.sugar==0] <- 0.1
red$chlorides[red$chlorides==0] <- 0.1
red$free.sulfur.dioxide[red$free.sulfur.dioxide==0] <- 0.1
red$total.sulfur.dioxide[red$total.sulfur.dioxide==0] <- 0.1
red$density[red$density==0] <- 0.1
red$pH[red$pH==0] <- 0.1
red$sulphates[red$sulphates==0] <- 0.1
red$alcohol[red$alcohol==0] <- 0.1
red$rspH[red$rspH==0] <- 0.1
red$catsd[red$catsd==0] <- 0.1
red$sa[red$sa==0] <- 0.1
red$fac[red$fac==0] <- 0.1
red$dpH[red$dpH==0] <- 0.1
red$fsda[red$fsda==0] <- 0.1

log_red <- log(red)

log_Red_stan_model_trans <- stan_glm(log(quality) ~ ., data = log_red, iter = 3000, warmup = 1500, prior_intercept = normal(1, 5), prior = normal(1, 1.5), adapt_delta = 0.999)
summary(log_Red_stan_model_trans)

mcmc_trace(as.array(log_Red_stan_model_trans), pars = c("(Intercept)", "volatile.acidity", 
                                                        "residual.sugar", "chlorides", "free.sulfur.dioxide" ,
                                                        "total.sulfur.dioxide", "pH", "sulphates", "rspH",
                                                        "catsd", "sa", "fac", "dpH", "fsda", "sigma"), n_warmup = 500)
#It doesn't converge


#Comparing models witH loo package
library(loo)

model_6 <- loo(Red_stan_model6)
model_select <- loo(Red_stan_model_select)
log_model <- loo(log_Red_stan_model)
model_trans <- loo(red_stan_model_trans)
model_trans_select <- loo(Red_stan_model_trans_sel)

loo_compare(model_6, model_select, log_model, model_trans, model_trans_select)











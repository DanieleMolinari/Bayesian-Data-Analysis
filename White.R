White <- read.csv("winequality-white.csv", header = TRUE, sep = ";")
White <- as.data.frame(White)
summary(White)
head(White)

hist(White$quality, xlab = "", main = "White wine quality", ylab = "Quality", breaks = 5, xlim = c(0, 10), col = "#4575b4")

correlWhite <- cor(White)
corrplot::corrplot(correlWhite, method = "number", type = "upper")

library(rstanarm)
library(bayesplot)

#Model1
set.seed(10)
White_stan_model <- stan_glm(quality ~ ., data = White, iter = 1000, warmup = 500)
summary(White_stan_model)

dimnames(as.array(White_stan_model))

color_scheme_set("blue")
mcmc_trace(as.array(White_stan_model), pars = c("(Intercept)", "sigma"),
           n_warmup = 500)

White_Rhats <- rhat(White_stan_model)
hat1 <- mcmc_rhat(White_Rhats)

#model 2
set.seed(10)
White_stan_model2 <- stan_glm(quality ~ ., data =White, iter = 1000, warmup = 500, prior_intercept = normal(0, 5), prior = normal(0, 2))
summary(White_stan_model2)

mcmc_trace(as.array(White_stan_model2), pars = c("(Intercept)", "sigma"),
           n_warmup = 500)

White_Rhats2 <- rhat(White_stan_model2)
hat2 <- mcmc_rhat(White_Rhats2)

#model 3
set.seed(10)
White_stan_model3 <- stan_glm(quality ~ ., data = White, iter = 1000, warmup = 500, prior_intercept = normal(0, 15), prior = normal(0, 5))
summary(White_stan_model3)

mcmc_trace(as.array(White_stan_model3), pars = c("(Intercept)", "sigma"),
           n_warmup = 500)

White_Rhats3 <- rhat(White_stan_model3)
hat3 <- mcmc_rhat(White_Rhats3)

#model 4
set.seed(10)
White_stan_model4 <- stan_glm(quality ~ ., data = White, iter = 1000, warmup = 500, prior_intercept = normal(1, 10), prior = normal(1, 2.5))
summary(White_stan_model4)


mcmc_trace(as.array(White_stan_model4), pars = c("(Intercept)", "sigma"),
           n_warmup = 500)

White_Rhats4 <- rhat(White_stan_model4)
hat4 <- mcmc_rhat(White_Rhats4)

#model 5
set.seed(10)
White_stan_model5 <- stan_glm(quality ~ ., data = White, iter = 1000, warmup = 500, prior_intercept = normal(1, 5), prior = normal(1, 2))
summary(White_stan_model5)

mcmc_trace(as.array(White_stan_model5), pars = c("(Intercept)", "sigma"),
           n_warmup = 500)

White_Rhats5 <- rhat(White_stan_model5)
hat5 <- mcmc_rhat(White_Rhats5)

#model 6
set.seed(10)
White_stan_model6 <- stan_glm(quality ~ ., data = White, iter = 1000, warmup = 500, prior_intercept = normal(1, 15), prior = normal(1, 5))
summary(White_stan_model6)

mcmc_trace(as.array(White_stan_model6), pars = c("(Intercept)", "fixed.acidity", "volatile.acidity", "citric.acid",
                                               "residual.sugar", "chlorides", "free.sulfur.dioxide" ,
                                               "total.sulfur.dioxide", "density", "pH", "sulphates", "alcohol", "sigma"), n_warmup = 500)

White_Rhats6 <- rhat(White_stan_model6)
hat6 <- mcmc_rhat(White_Rhats6)

library(cowplot)
plot_grid(hat1, hat2, hat3, hat4, hat5, hat6, ncol = 3, labels = c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6"), hjust = -2)

plot(White_stan_model5, "acf", pars = c("(Intercept)", "fixed.acidity", "volatile.acidity", "citric.acid",
                                      "residual.sugar", "chlorides", "free.sulfur.dioxide" ,
                                      "total.sulfur.dioxide", "density", "pH", "sulphates", "alcohol", "sigma"))
mcmc_trace(as.array(White_stan_model5), pars = c("(Intercept)", "fixed.acidity", "volatile.acidity", "citric.acid",
                                                 "residual.sugar", "chlorides", "free.sulfur.dioxide" ,
                                                 "total.sulfur.dioxide", "density", "pH", "sulphates", "alcohol", "sigma"), n_warmup = 500)


#############################################################################
################Lasso with Metropolis - Hastings sampler#####################

XWhite <- White[,-12]
yWhite <- White$quality
XWhite <- matrix(unlist(XWhite), ncol = dim(White)[2], byrow = FALSE)
XWhite <- cbind(rep(1, nrow(XWhite)), scale(XWhite))
n <- length(yWhite)
k <- ncol(XWhite)
npars <- k + (k-1) + 2
nchains <- 3

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
             dgamma(1/sigma2,1,15,log=TRUE) +
             dgamma(lambda2,0.02,0.1,log=TRUE)+
             sum(dexp(tau,lambda2/2,log=TRUE)))
  return(dens)
}

loglikelihood <- function(yWhite,XWhite,par) {
  b <- as.matrix(par[1:k])
  sigma2 <- par[npars]
  if (sigma2 <= 0)
    return(log(0))
  dens <- sum(dnorm(yWhite,t(XWhite%*%b),sigma2,log=TRUE))
  return(dens)
}

draws <- matrix(0,nrow=40000,ncol=npars)
draws[1,npars] <- 1/rgamma(1,1,15)
draws[1,npars-1] <- 1
draws[1,1] <- mean(yWhite)
draws[1,2:k] <- rep(0,k-1)
draws[1,(k+1):(2*k-1)] <- rep(1,k-1)
for (step in 2:40000) {
  proposed <- rnorm(npars,draws[step-1,],0.1)
  r <- loglikelihood(yWhite,XWhite,proposed)+
    logprior(proposed)-
    loglikelihood(yWhite,XWhite,draws[step-1,])-
    logprior(draws[step-1,])
  u <- runif(1)
  if (log(u) < r) {
    draws[step,] <- proposed
  } else {
    draws[step,] <- draws[step-1,]
  }
}
sample <- draws[20000:40000,]
betas <- sample[, 3:k]

boxplot(betas, outline = FALSE, horizontal = TRUE)
abline(v = 0, lty = 2)
#############################################################################
#############################################################################


#Blasso
library(monomvn)
yWhite <- White$quality
XWhite <- White[,-12]
XWhite <- matrix(unlist(XWhite), ncol = dim(XWhite)[2], byrow = FALSE)
XWhite <- scale(XWhite)

Whitebla <- blasso(XWhite, yWhite, T=11000, rd = c(1, 5), ab = c(1, 2),
                   thin = 100, RJ = FALSE, normalize = FALSE)

White_plot <- boxplot(Whitebla$beta, outline = FALSE, horizontal = TRUE)
abline(v = 0, lty = 2)
#take away 3, 5, 7

#select model
White_fin <- White[, -c(3, 5, 7)]
set.seed(10)
White_stan_model_select <- stan_glm(quality ~ ., data = White_fin, iter = 1000, warmup = 500, prior_intercept = normal(1, 5), prior = normal(1, 2))
mcmc_trace(as.array(White_stan_model_select), pars = c("(Intercept)","fixed.acidity", "volatile.acidity", "residual.sugar", "free.sulfur.dioxide",
                                                       "density", "pH", "sulphates", "alcohol" , "sigma"), n_warmup = 500)

White_Rhats_select <- rhat(White_stan_model_select)
hat_select <- mcmc_rhat(White_Rhats_select)
hat_select

plot(White_stan_model_select, "acf", pars = c("(Intercept)","fixed.acidity", "volatile.acidity", "residual.sugar", "free.sulfur.dioxide",
                                              "density", "pH", "sulphates", "alcohol" , "sigma"))

#predictions with model 5
predictions_firstmodel <- posterior_linpred(White_stan_model5)

iteration_1_5<- predictions_firstmodel[1,]
iteration_328_5 <- predictions_firstmodel[328,]
iteration_500_5 <- predictions_firstmodel[500,]
iteration_1025_5 <- predictions_firstmodel[1025,]


summary(White$quality)

summary(iteration_1_5)
summary(iteration_328_5)
summary(iteration_500_5)
summary(iteration_1025_5)

yrep_5 <- posterior_predict(White_stan_model5)

par(mfrow = c(2, 3))

first_5 <- predictions_firstmodel[,252]
summary(first_5)
White$quality[252]
plot(density(yrep_6[,252]), col = "red")
abline(v = White[252, 12])

second_5 <- predictions_firstmodel[,2819]
summary(second_5)
White[2819,12]
plot(density(yrep_6[,2819]), col = "red")
abline(v = White[2819, 12])

third_5 <- predictions_firstmodel[,268]
summary(third_5)
White[268,12]
plot(density(yrep_6[,268]), col = "red")
abline(v = White[268, 12])

fourth_5 <- predictions_firstmodel[,487]
summary(fourth_5)
White[487,12]
plot(density(yrep_6[,487]), col = "red")
abline(v = White[487, 12])

fifth_5 <- predictions_firstmodel[, 1163]
summary(fifth_5)
White[1163,12]
plot(density(yrep_6[,1163]), col = "red")
abline(v = White[1163, 12])

sixth_5 <- predictions_firstmodel[,3030]
summary(sixth_5)
White[3030,12]
plot(density(yrep_6[,3030]), col = "red")
abline(v = White[3030, 12])

pp_3 <- pp_check(White_stan_model3, "dens_overlay", adjust = 5)
pp_check(White_stan_model3, "stat")
pp_check(White_stan_model3, "stat_2d")

#Prediction model selected
predictions_selectmodel <- posterior_linpred(White_stan_model_select)

yrep_sel <- posterior_predict(White_stan_model_select)

first_Sel <- predictions_selectmodel[,252]
summary(first_Sel)
White$quality[252]
plot(density(yrep_sel[,252]), col = "red")
abline(v = White[252, 12])

second_Sel <- predictions_selectmodel[,2819]
summary(second_Sel)
White[2819,12]
plot(density(yrep_sel[,2819]), col = "red")
abline(v = White[2819, 12])

third_Sel <- predictions_selectmodel[,268]
summary(third_Sel)
White[268,12]
plot(density(yrep_sel[,268]), col = "red")
abline(v = White[268, 12])

fourth_Sel <- predictions_selectmodel[,487]
summary(fourth_Sel)
White[487,12]
plot(density(yrep_sel[,487]), col = "red")
abline(v = White[487, 12])

fifth_Sel <- predictions_selectmodel[,1163]
summary(fifth_Sel)
White[1163,12]
plot(density(yrep_sel[,1163]), col = "red")
abline(v = White[1163, 12])

sixth_Sel <- predictions_selectmodel[,3030]
summary(sixth_Sel)
White[3030 ,12]
plot(density(yrep_sel[, 3030]), col = "red")
abline(v = White[3030, 12])

pp_sel1 <- pp_check(White_stan_model_select, "dens_overlay", adjust = 5)
pp_check(White_stan_model_select, "stat")
pp_check(White_stan_model_select, "stat_2d")

library(cowplot)
plot_grid(pp_3, pp_sel1, labels = c("Model 3", "Model Selected"),  hjust = -1.75)

#Do the log transormation
White$fixed.acidity[White$fixed.acidity==0] <- 0.1
White$volatile.acidity[White$volatile.acidity==0] <- 0.1
White$citric.acid[White$citric.acid==0] <- 0.1
White$residual.sugar[White$residual.sugar==0] <- 0.1
White$chlorides[White$chlorides==0] <- 0.1
White$free.sulfur.dioxide[White$free.sulfur.dioxide==0] <- 0.1
White$total.sulfur.dioxide[White$total.sulfur.dioxide==0] <- 0.1
White$density[White$density==0] <- 0.1
White$pH[White$pH==0] <- 0.1
White$sulphates[White$sulphates==0] <- 0.1
White$alcohol[White$alcohol==0] <- 0.1

log_White <- log(White)

log_White_stan_model <- stan_glm(log(quality) ~ log(fixed.acidity) + log(volatile.acidity) + log(citric.acid) + 
                                   log(residual.sugar) + log(chlorides) + log(free.sulfur.dioxide) + 
                                   log(total.sulfur.dioxide) + log(density) + log(pH) + log(sulphates) +
                                   log(alcohol), data = White, iter = 1000, warmup = 500, prior_intercept = normal(1, 5), prior = normal(1, 2), adapt_delta = 0.999)
summary(log_White_stan_model)

mcmc_trace(as.array(log_White_stan_model), pars = c("(Intercept)", "fixed.acidity", "volatile.acidity", "citric.acid",
                                                  "residual.sugar", "chlorides", "free.sulfur.dioxide" ,
                                                  "total.sulfur.dioxide", "density", "pH", "sulphates", "alcohol", "sigma"), n_warmup = 500)

White_Rhats_log_model <- rhat(log_White_stan_model)
hat_log_model <- mcmc_rhat(White_Rhats_log_model)
hat_log_model

plot(log_White_stan_model, "acf", pars = c("(Intercept)", "fixed.acidity", "volatile.acidity", "citric.acid",
                                           "residual.sugar", "chlorides", "free.sulfur.dioxide" ,
                                           "total.sulfur.dioxide", "density", "pH", "sulphates", "alcohol", "sigma"))

#predictions with white log model
predictions_log_model <- posterior_linpred(log_White_stan_model)

iteration_1_log <- predictions_log_model[1,]
iteration_328_log <- predictions_log_model[328,]
iteration_500_log <- predictions_log_model[500,]
iteration_1025_log <- predictions_log_model[1025,]

summary(log(White$quality))

summary(iteration_1_log)
summary(iteration_328_log)
summary(iteration_500_log)
summary(iteration_1025_log)

yrep_log <- posterior_predict(log_White_stan_model)

par(mfrow = c(2, 3))

first_log <- predictions_log_model[,252]
summary(first_log)
log(White$quality[252])
plot(density(yrep_log[,252]), col = "red",  main = "Quality Wine 3")
abline(v = log(White[252, 12]))

second_log <- predictions_log_model[,2819]
summary(second_log)
log(White[2819,12])
plot(density(yrep_log[,2819]), col = "red",  main = "Quality Wine 4")
abline(v = log(White[2819, 12]))

third_log <- predictions_log_model[,268]
summary(third_log)
log(White[268, 12])
plot(density(yrep_log[, 268]), col = "red",  main = "Quality Wine 5")
abline(v = log(White[268, 12]))

fourth_log <- predictions_log_model[,487]
summary(fourth_log)
log(White[487,12])
plot(density(yrep_log[, 487]), col = "red",  main = "Quality Wine 6")
abline(v = log(White[487, 12]))

fifth_log <- predictions_log_model[,1163]
summary(fifth_log)
log(White[1163, 12])
plot(density(yrep_log[, 1163]), col = "red",  main = "Quality Wine 7")
abline(v = log(White[1163, 12]))

sixth_log <- predictions_log_model[,3030]
summary(sixth_log)
log(White[3030, 12])
plot(density(yrep_log[, 3030]), col = "red",  main = "Quality Wine 8")
abline(v = log(White[1270, 12]))

pp_check(log_White_stan_model, "dens_overlay", adjust = 5)
pp_check(log_White_stan_model, "stat")
pp_check(log_White_stan_model, "stat_2d")


#Add pairwise variables
White$rspH <- White$residual.sugar*White$pH
White$catsd <- White$citric.acid*White$total.sulfur.dioxide
White$sa <- White$sulphates*White$alcohol
White$fac <- White$fixed.acidity/White$chlorides
White$dpH <- White$density/White$pH
White$fsda <- White$free.sulfur.dioxide/White$alcohol

set.seed(10)
White_stan_model_trans <- stan_glm(quality ~ ., data = White, iter = 1000, warmup = 500, prior_intercept = normal(1, 5), prior = normal(1, 2))
summary(White_stan_model_trans)

mcmc_trace(as.array(White_stan_model_trans), pars = c("(Intercept)", "fixed.acidity", "volatile.acidity", "citric.acid",
                                                      "residual.sugar", "chlorides", "free.sulfur.dioxide" ,
                                                      "total.sulfur.dioxide", "density", "pH", "sulphates", "alcohol", "rspH",
                                                      "catsd", "sa", "fac", "dpH", "fsda", "sigma"), n_warmup = 500)

White_Rhats_model_trans <- rhat(White_stan_model_trans)

mcmc_rhat(White_Rhats_model_trans)

predictions_trans <- posterior_linpred(White_stan_model_trans)

yrep_trans <- posterior_predict(White_stan_model_trans)

first_trans <- predictions_trans[,252]
summary(first_trans)
White$quality[252]
plot(density(yrep_trans[,252]), col = "red")
abline(v = White[252, 12])

second_trans <- predictions_trans[,2819]
summary(second_trans)
White[2819,12]
plot(density(yrep_trans[,2819]), col = "red")
abline(v = White[2819, 12])

third_trans <- predictions_trans[, 268]
summary(third_trans)
White[268, 12]
plot(density(yrep_trans[, 268]), col = "red")
abline(v = White[268, 12])

fourth_trans <- predictions_trans[,487]
summary(fourth_trans)
White[487,12]
plot(density(yrep_trans[, 487]), col = "red")
abline(v = White[487, 12])

fifth_trans <- predictions_trans[, 1163]
summary(fifth_trans)
White[1163, 12]
plot(density(yrep_trans[, 1163]), col = "red")
abline(v = White[1163, 12])

sixth_trans <- predictions_trans[, 3030]
summary(sixth_trans)
White[3030, 12]
plot(density(yrep_trans[, 3030]), col = "red")
abline(v = White[3030, 12])

par(mfrow = c(1, 1))

pp_check(White_stan_model_trans, "dens_overlay", adjust = 5)
pp_check(White_stan_model_trans, "stat")
pp_check(White_stan_model_trans, "stat_2d")


#Reapplying lasso scratch...

#Reapplying Blasso function

yWhite <- White$quality
XWhite <- White[,-12]
XWhite <- matrix(unlist(XWhite), ncol = dim(XWhite)[2], byrow = FALSE)
XWhite <- scale(XWhite)

Whitebla <- blasso(XWhite, yWhite, T=11000, rd = c(1, 5), ab = c(0, 2),
                   thin = 100, RJ = FALSE, normalize = FALSE)

White_plot <- boxplot(Whitebla$beta, outline = FALSE, horizontal = TRUE)
abline(v = 0, lty = 2)

#Taking away 3, 5, 7, 9, 10, 16

White_fin_2 <- White[, -c(3, 5, 7, 9, 10, 16)]
set.seed(10)
White_stan_model_trans_sel <- stan_glm(quality ~ ., data = White_fin_2, iter = 1000, warmup = 500, prior_intercept = normal(1, 5), prior = normal(1, 2))

mcmc_trace(as.array(White_stan_model_trans_sel), pars = c("(Intercept)","fixed.acidity", "volatile.acidity", "residual.sugar", "free.sulfur.dioxide",
                                                       "density", "alcohol","rspH", "catsd", "sa", "dpH", "fsda" , "sigma"), n_warmup = 500)


White_Rhats_model_trans_sel <- rhat(White_stan_model_trans_sel)

mcmc_rhat(White_Rhats_model_trans_sel)

predictions_trans_sel <- posterior_linpred(White_stan_model_trans_sel)

yrep_trans_sel <- posterior_predict(White_stan_model_trans_sel)

first_trans_sel <- predictions_trans_sel[,252]
summary(first_trans_sel)
White$quality[252]
plot(density(yrep_trans_sel[,252]), col = "red")
abline(v = White[252, 12])

second_trans_sel <- predictions_trans_sel[,2819]
summary(second_trans_sel)
White[2819,12]
plot(density(yrep_trans_sel[,2819]), col = "red")
abline(v = White[2819, 12])

third_trans_sel <- predictions_trans_sel[, 268]
summary(third_trans_sel)
White[268, 12]
plot(density(yrep_trans_sel[, 268]), col = "red")
abline(v = White[268, 12])

fourth_trans_sel <- predictions_trans_sel[,487]
summary(fourth_trans_sel)
White[487,12]
plot(density(yrep_trans_sel[, 487]), col = "red")
abline(v = White[487, 12])

fifth_trans_sel <- predictions_trans_sel[, 1163]
summary(fifth_trans_sel)
White[1163, 12]
plot(density(yrep_trans_sel[, 1163]), col = "red")
abline(v = White[1163, 12])

sixth_trans_sel <- predictions_trans_sel[, 3030]
summary(sixth_trans_sel)
White[3030, 12]
plot(density(yrep_trans_sel[, 3030]), col = "red")
abline(v = White[3030, 12])

pp_check(White_stan_model_trans_sel, "dens_overlay", adjust = 5)
pp_check(White_stan_model_trans_sel, "stat")
pp_check(White_stan_model_trans_sel, "stat_2d")


#Comparing models with package loo

model_5 <- loo(White_stan_model3)
model_select_w <- loo(White_stan_model_select)
log_model_w <- loo(log_White_stan_model)
model_trans_w <- loo(White_stan_model_trans)
model_trans_select_w <- loo(White_stan_model_trans_sel)

loo_compare(model_5, model_select_w, log_model_w, model_trans_w, model_trans_select_w)






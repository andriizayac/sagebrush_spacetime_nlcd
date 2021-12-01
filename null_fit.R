pkgs <- c("dplyr", "tidyr", "lme4", "ggplot2", "brms")
sapply(pkgs, require, character.only = T)

years <- c(1985:2018)

source("helper_fns.R")
# ====================================
# see pxlmatching.R 
tfires <- readRDS("data/tfires.rds")
tsage <- readRDS("data/tsage.rds")
tpxlcov <- readRDS("data/tpxlcov.rds")
tdfEnv <- readRDS("data/tdfEnv_covars.rds")
txys <- readRDS("data/txys.rds")

N <- nrow(tfires)

kvec <- sapply(tpxlcov, function(x){ mean(x$prefire)})

# --- select 5% of the observations from each sample 
ssize <- sapply(tsage, function(x) floor(nrow(x)*.1) )
nullsage <- lapply(tsage, function(x) {
  ind = sample(1:nrow(x), floor(nrow(x)*.1), replace = F)
  x[ind,]
})

# --- combine data sets - different number of years post-fire don't matter because model is fit recursively
dat.gen.null.full <- function(datin, tfires, i) {
  mat <- nullsage[[i]][,c(tfires$FireYer[i]-1984):31]
  T <- ncol(mat)
  colnames(mat) <- 1:T-1
  dat <- mat %>% 
    pivot_longer(cols = 1:T) %>%
    rename(t = name, cover = value) %>% 
    mutate(t = as.numeric(t)) 
  df <- data.frame(y = dat$cover[dat$t != 0],
                   x = dat$cover[dat$t != max(dat$t)],
                   t = dat$t[dat$t != max(dat$t)], 
                   id = i)
  return(df)
}
dat.gen.null <- function(datin, tfires, i) {
  mat <- nullsage[[i]][,c(tfires$FireYer[i]-1984):31]
  T <- ncol(mat)
  colnames(mat) <- 1:T-1
  dat <- mat %>% 
    pivot_longer(cols = 1:T) %>%
    rename(t = name, cover = value) %>% 
    mutate(t = as.numeric(t)) 
  df <- data.frame(y = dat$cover[dat$t != 0],
                   x = dat$cover[dat$t != max(dat$t)],
                   t = dat$t[dat$t != max(dat$t)], 
                   id = i) %>% 
    filter(x != 0)
    return(df)
}
dat.gen.null.init <- function(datin, tfires, i=NULL){
  mat <-datin[[i]][,c(tfires$FireYer[i]-1983):31]
  T = ncol(mat)
  colnames(mat) <- 1:T-1
  dat <- pivot_longer(mat, cols = c(1:T)) %>%
    rename(t = name, cover = value) %>% 
    mutate(t = as.numeric(t)) %>% 
    filter(t <= 5) %>% 
    mutate(id = i)
  return(dat)
}

# --- initiate a combined data set and append the rest
dff <- dat.gen.null.full(nullsage, tfires, 1)
df <- dat.gen.null(nullsage, tfires, 1)
df0 <- dat.gen.null.init(nullsage, tfires, 1)
for(i in 2:N) {
  dff <- rbind(dff, dat.gen.null.full(nullsage, tfires, i))
  df <- rbind(df, dat.gen.null(nullsage, tfires, i))
  df0 <- rbind(df0, dat.gen.null.init(nullsage, tfires, i))
}

# --- estimate growth and dd parameters using lme4
temp <- glm(y ~ 1 + log(x), offset = log(x), family = poisson(link = "log"), data = df) #
# --- estimate initial population size based the first 5 years of data post-fire
temp_n0 <- glm(cover ~ 1 + t, family = poisson(link = "log"), data = df0) # 

m <- coef(temp)
k <- exp(-m[1]/m[2])
n0 <- exp(coef(temp_n0)[1])

t <- seq(0, 33, by = .1)
pred <- gomp(m[1], m[2], n0, t)[1,] 
plot(pred ~ t, type = "l", lwd = 4, ylim = c(0, 25),
     main = "Average predictions across the GB", 
     ylab = "Predicted cover, %", xlab = "Time, years")
abline(h = k, lty = "dashed", col = "gray", lwd = 1.5)
# --- add real data
lines(r$y ~ r$t, lwd = 4, lty = "dashed", col = "brown")
#apply(predb[3000:3100,], 1, function(x) {lines(x, type = "l", col = rgb(0,0,0,.1))})
cat("Region-wide 'equilibrium' abundance is:", k, "%")

# === MAE under the Gompertz model in t using lme4 modles
ymus <- list()
yhats <- list()

maenull <- rep(NA, N)
for(i in 1:N) {
  y <- nullsage[[i]][, c(tfires$FireYer[i]-1984):31] %>% 
    summarize_all(mean) %>% 
    as.matrix()
  d <- dim(y)
  yhat <- matrix(gomp(m[1], m[2], n0, 1:d[2]), 
                 nr = d[1], 
                 nc = d[2], byrow = T)
  ymus[[i]] <- y
  yhats[[i]] <- yhat
  # maenull[i] <- mean(abs(y - yhat))
}
plot(pred ~ t, type = "l", lwd = 4, ylim = c(0, 25),
     main = "Average predictions across the GB", 
     ylab = "Predicted cover, %", xlab = "Time, years")
sapply(ymus, function(x) { lines(1:length(x), x, lwd = 1, col = rgb(0,1,1,.2)) } )
lines(pred ~ t, type = "l", lwd = 4)
abline(h = k, lty = "dashed", col = "gray", lwd = 1.5)

# plot(density(maenull), lwd = 2, main = "MAE Null model", xlab = "%, abundance", xlim = c(0,20))
# === MAE using the recursive form Gompertz (in-sample only!)
# --- initiate with the first case
# mat <- data.matrix(nullsage[[1]][,c(tfires$FireYer[1]-1984):33])
# xpr <- vec(mat[, -ncol(mat)])
# ypr <- vec(mat[, -1])
# yhat <- predict.glm(temp, newdata = exp(data.frame(x=xpr)))
# 
# inrecerr <- rep(NA, N)
# inrecerr[1] <- mean(abs(yhat - ypr))
# for(i in 2:N) {
#   svMisc::progress(i)
#   mat <- data.matrix(nullsage[[i]][,c(tfires$FireYer[i]-1984):33])
#   xpr <- c(xpr, vec(mat[, -ncol(mat)]))
#   ypr <- c(ypr, vec(mat[, -1]))
#   yhat <- predict.glm(temp, newdata = exp(data.frame(x=xpr)))
#   inrecerr[i] <- mean(abs(yhat - ypr))
# }
# 
# lines(density(inrecerr), lwd = 2, lty = "dashed", col = "lightblue", main = "MAE Null - in sample", xlab = "%, abundance")
# 
# data.frame(recursive = inrecerr, gomp = mae) %>% 
#   pivot_longer(cols = 1:2) %>% 
#   ggplot(aes(x = name, y = value)) + geom_boxplot() +
#   labs(x = "Model", y = "MAE, %") +
#   theme_bw()
# =====


plot(maenull ~kvec, pch = 19, col = rgb(0,0,0,.5),
     ylab = "MEA null", xlab = "Mean pre-disturbance cover, %")
abline(0, 1)

# saveRDS(list(temp = temp, temp_n0 = temp_n0, 
#              df = df, df0 = df0, dffull = dff, 
#              maenull = maenull, ymus = ymus), 
#         file = "outputs/null_model.rds")
# nulldat <- readRDS("outputs/null_model.rds")

# === brms
######################################################
nulldat <- readRDS("outputs/null_dat.rds")
ymus <- list()
for(i in 1:N) {
  ymus[[i]] <- nulldat[[i]][, c(tfires$FireYer[i]-1984):31] %>% 
    summarize_all(mean) %>% 
    as.matrix()
}
# === using Bayesian estimation to get pars uncertainty
library(brms)
# input <- readRDS("outputs/null_model.rds")
# bdat <- input$df
# bdat0 <- input$df0
# --- fit brms models
# tempb <- brm(y ~ 1 + log(x) + offset(log(x)),
#            data = bdat,
#            family = poisson, 
#            algorithm = "meanfield")
# 
# tempb0 <- brm(cover ~ 1 + t,
#              data = bdat0,
#              family = poisson, 
#              algorithm = "meanfield")
# ---
tempb <- readRDS("outputs/null_model_brm.rds")
tempb0 <- readRDS("outputs/null_model0_brm.rds")
# ---
cb <- posterior_samples(tempb)
cb0 <- exp( posterior_samples(tempb0)$b_Intercept )
mb <- exp(-cb$b_Intercept / cb$b_logx)
t <- 0:50

predb <- gomp(cb$b_Intercept, cb$b_logx, cb0, t)
predb.mu <- apply(predb, 2, mean)
predb.ci <- apply(predb, 2, quantile, probs = c(0.025, 0.975))

# --- mae
maenull <- matrix(NA, 4000, N)
for(i in 1:N) {
  y <- ymus[[i]][1,9]
  maenull[, i] = sapply(predb[,9], function(x) {mean(abs(y - x))} )
}

# ---
# pdf("figures/fig2.pdf", width = 7, height = 3.5)
par(mfrow = c(1, 2), mar = c(4.5,4.5,1,1))
plot(1, type = "n", xlim = range(t), ylim = c(0, 18), ylab = "Predicted cover, %", xlab = "Time, years")
sapply(ymus, function(x) { lines(1:length(x)-1, x, lwd = 1, col = rgb(0,1,1,.2)) } )
# polygon(c(t, rev(t)), c(predb.ci[2, ], rev(predb.ci[1, ])), col = rgb(.5,.5,.5,.8), border = FALSE)
apply(predb, 1, function(x) {lines(x ~ t, col = rgb(.6,0.6,0.6,1)) } )
lines(t, predb.mu, lty = "dashed", lwd = 4, col = "black")

polygon(c(0, max(t), max(t), 0), sort(rep(quantile(mb, probs = c(0.025, 0.975)), 2)), col = rgb(.5,.5,.5,.8), border = FALSE)
lines(t, rep(mean(mb), length(t)), lty = "dashed", lwd = 2, col = "black")

hist(maenull, main = "", xlab = "MAE, %", col = rgb(.6,0.6,0.6,1))
# dev.off()

mae.null.mu <- apply(maenull, 2, mean)
mae.null.ci <- apply(maenull, 2, quantile, prob = c(0.025, 0.975))
plot(density(mae.null.mu), lwd = 3, main = "", xlab = "MAE, %")
lines(density(mae.null.ci[1, ]), lwd = 2, col = "red")
lines(density(mae.null.ci[2, ]), lwd = 2, col = "red")
# === end of bayes estimate



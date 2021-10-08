pkgs <- c("dplyr", "tidyr", "lme4", "ggplot2")
sapply(pkgs, require, character.only = T)

years <- c(1985:2018)

source("helper_fns.R")
# ====================================
# see pxlmatching.R 
tfires <- readRDS("data/tfires.rds")
tsage <- readRDS("data/tsage.rds")
tpxlcov <- readRDS("data/tpxlcov.rds")
tdfEnv <- readRDS("data/tdfEnv_covars.rds")

N <- nrow(tfires)

# --- select 5% of the observations from each sample 
ssize <- sapply(tsage, function(x) floor(nrow(x)*.05) )
nullsage <- lapply(tsage, function(x) {
  ind = sample(1:nrow(x), floor(nrow(x)*.05), replace = F)
  x[ind,]
})

# --- combine data sets - different number of years post-fire don't matter because model is fit recursively
dat.gen.null <- function(datin, tfires, i) {
  mat <- data.matrix(nullsage[[i]][,c(tfires$FireYer[i]-1984):33])
  xpr <- vec(mat[, -ncol(mat)])
  ypr <- vec(mat[, -1])
  outs <- which(xpr == 0 | ypr == 0)
  dat <- data.frame(y = ypr[-outs],
                    x = xpr[-outs]) 
}
dat.gen.null.init <- function(tsage, tfires, i=NULL){
  mat <- tsage[[i]][,c(tfires$FireYer[i]-1984):33]
  T = ncol(mat)
  colnames(mat) <- 1:T-1
  dat <- pivot_longer(mat, cols = c(1:T)) %>%
    rename(t = name, cover = value) %>% 
    mutate(t = as.numeric(t), 
           cover = ifelse(cover == 0, 0, cover)) %>% #runif(1e6, 0,1)
    filter(t <= 5)
  return(dat)
}
# --- initiate a combined data set and append the rest
df <- dat.gen.null(nullsage, tfires, 1)
df0 <- dat.gen.null.init(nullsage, tfires, 1)
for(i in 2:N) {
  df <- rbind(df, dat.gen.null(nullsage, tfires, i))
  df0 <- rbind(df0, dat.gen.null.init(nullsage, tfires, i))
}

# --- estimate growth and DD parameters
temp <- glm(y ~ 1 + log(x), offset = log(x), family = poisson(), data = df) #
# --- estimate initial population size based the first 5 years of data post-fire
temp_n0 <- glm(cover ~ 1 + t, family = poisson(link = "log"), data = df0) # 
# === using Bayesian estimates to estimate uncertainty
# library(brms)
# set.seed(124)
# idx <- sample(1:nrow(df), 10000)
# tempb <- brm(y ~ 1 + log(x) + offset(log(x)),  
#            data = df,#[idx,],
#            family = poisson)
# idx0 <- sample(1:nrow(df0), 10000)
# tempb0 <- brm(cover ~ 1 + t,  
#              data = df0,#[idx0,],
#              family = poisson)
# cb <- posterior_samples(tempb)
# cb0 <- posterior_samples(tempb0)
# mb <- exp(-cb[,1]/cb[,2])
# predb <- gomp(cb[,1], cb[,2], exp(cb[,1]), t)
# predb <- apply(predb, 2, quantile, probs = c(0.025, 0.975))
# === end of bayes estimate

m <- coef(temp)
k <- exp(-m[1]/m[2])
n0 <- exp(coef(temp_n0)[1])

t <- seq(1, 70, by = .1)
pred <- gomp(m[1], m[2], n0, t)[1,] 
plot(pred ~ t, type = "l", lwd = 2, 
     main = "Average predictions across the GB", 
     ylab = "Predicted cover, %", xlab = "time, years")
abline(h = k, lty = "dashed", col = "gray", lwd = 1.5)
#apply(predb[3000:3100,], 1, function(x) {lines(x, type = "l", col = rgb(0,0,0,.1))})
cat("Region-wide 'equilibrium' abundance is:", k, "%")
# === MAE using the Gompertz model in t
mae <- rep(NA, N)
for(i in 1:N) {
  y <- tsage[[i]] %>% as.matrix()
  d <- dim(y)
  yhat <- matrix(gomp(m[1], m[2], n0, 1:d[2]), 
                 nr = d[1], 
                 nc = d[2], byrow = T)
  mae[i] <- mean(abs(y - yhat))
}
plot(density(mae), lwd = 2, main = "MAE Null model", xlab = "%, abundance", xlim = c(0,20))
# === MAE using the recursive form Gompertz (in-sample only!)
# --- initiate with the first case
mat <- data.matrix(nullsage[[1]][,c(tfires$FireYer[1]-1984):33])
xpr <- vec(mat[, -ncol(mat)])
ypr <- vec(mat[, -1])
yhat <- predict.glm(temp, newdata = exp(data.frame(x=xpr)))

inrecerr <- rep(NA, N)
inrecerr[1] <- mean(abs(yhat - ypr))
for(i in 2:N) {
  svMisc::progress(i)
  mat <- data.matrix(nullsage[[i]][,c(tfires$FireYer[i]-1984):33])
  xpr <- c(xpr, vec(mat[, -ncol(mat)]))
  ypr <- c(ypr, vec(mat[, -1]))
  yhat <- predict.glm(temp, newdata = exp(data.frame(x=xpr)))
  inrecerr[i] <- mean(abs(yhat - ypr))
}

lines(density(inrecerr), lwd = 2, lty = "dashed", col = "lightblue", main = "MAE Null - in sample", xlab = "%, abundance")

data.frame(recursive = inrecerr, gomp = mae) %>% 
  pivot_longer(cols = 1:2) %>% 
  ggplot(aes(x = name, y = value)) + geom_boxplot() +
  labs(x = "Model", y = "MAE, %") +
  theme_bw()



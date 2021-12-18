pkgs <- c("dplyr", "tidyr")
sapply(pkgs, require, character.only = T)

years <- c(1985:2018)

source("helper_fns.R")
tpxlcov <- readRDS("data/tpxlcov.rds")
tdfEnv <- readRDS("data/tdfEnv_covars.rds")
tfires <- readRDS("data/tfires.rds")
tsage <- readRDS("data/tsage.rds")


coefs_brm <- t( sapply(louts, function(x)
                    apply(x, 2, mean)) )
# --- site level
# modout <- readRDS("nlcd/model_output_test/modout_1.rds")
# datlist <- readRDS("nlcd/model_output_test/datlist_1.rds")
modout <- readRDS("~/Downloads/modout_1.rds")
datlist <- readRDS("~/Downloads/datlist_1.rds")

coefs <- modout$coefs
tvec <- modout$tvec
  
N <- length(coefs)

# tdfEnv$ninit <- sapply(coefs, function(x) { x$n0 } )
#wt <- diag( c( 1, 1, 1, 1, 0.4283099, 0.4494240, 0.4315081, 3.4884838, 3.0928078, 0.8467186 ) )
wt <- diag( c( 1, 1, 1, 1, 1.745238,  6.673595,  8.504045, 13.378898,  8.391983,  8.406849 ) )
# wt <- diag( rep(1, 10))
cov <- var(tdfEnv) %*% wt
mah <- StatMatch::mahalanobis.dist(tdfEnv, vc = cov)
diag(mah) <- NA
outs <- which( sapply(coefs, function(x) { any(x$b > 0) } ) == TRUE)
# outs2 <- which( sapply(coefs, function(x) { mean(x$n0) > 5 } ) == TRUE)
mah[, outs] <- NA
# mah[, outs2] <- NA
match <- apply(mah, 1, which.min)

matchsite <- list(pred = list(), obs = list())

for(i in 1:N){
  print(i)
  obs <- apply(datlist[[i]], 2, mean)
  # out-of-sample Mahalanobis nearest match
  j <- match[i]
  a <- coefs[[j]]$a
  b <- coefs[[j]]$b
  n0 <- coefs[[j]]$n0
  t <- 1:tvec[i]-1
  matpred.m <- matrix(gomp(a, b, n0, t), 
                      nr = 1, #nrow(obs), 
                      nc = length(t), byrow = T)
  
  matchsite$pred[[i]] <- matpred.m
  matchsite$obs[[i]] <- obs
}

#saveRDS(matchsite, file = paste0("nlcd/model_pred/match_pred_1.rds"))
m <- rep(NA, N)
for(i in 1:N){
  yhat <- matchsite$pred[[i]]
  m[i] = mean(abs( matchsite$obs[[i]][1:length(yhat)] - yhat ))
}

# --- site + cluster levels
lout <- list()
clvec <- 2:14
for(M in c(4, 14)){
cat("---working on -", M, "-cluster---")  
  # inm <- paste0("nlcd/model_output_test/modout_", M, ".rds")
  # ind <- paste0("nlcd/model_output_test/datlist_", M, ".rds")
  inm <- paste0("~/Downloads/modout_", M, ".rds")
  ind <- paste0("~/Downloads/datlist_", M, ".rds")
  
  modout <- readRDS(inm)
  datlist <- readRDS(ind)
  
  coefs <- modout$coefs
  tvec <- modout$tvec
  
  msing <- unlist(modout$msing)
  msing0 <- unlist(modout$msing0)
  
  N <- length(coefs)
# === InSample predictions
##################################
#insample <- list(pred = list())
#names(errlist) <- paste0("cluster", 3:15)
#for(j in 1:N) {
#	a <- coefs[[j]]$a
#	b <- coefs[[j]]$b
#	n0 <- coefs[[j]]$n0
#	t <- 1:tvec[j]-1
#	matpred <- gomp(a, b, n0, t) %>%
#		as.data.frame() %>%
#		mutate(cl = row_number())
#
#	jmat <- datlist[[j]] %>%
#		as.data.frame() %>%
#		mutate(cluster = tpxlcov[[j]][, paste0("cluster", M)]) %>%
#		left_join(matpred, by = c("cluster" = "cl"))
#	pred <- jmat %>% select(starts_with("V")) %>% as.matrix()
#	insample$pred[[j]] <- pred
#}
#saveRDS(insample, file = paste0("nlcd/model_pred/insample_pred_", M, ".rds"))

# === Matching at site and pixel levels (pixel-to-cluster)
# tdfEnv$ninit <- sapply(coefs, function(x) { x$n0 } )
wt <- diag( c( 1, 1, 1, 1, 0.4283099, 0.4494240, 0.4315081, 3.4884838, 3.0928078, 0.8467186 ) )
# wt <- diag( rep(1, 10))
cov <- var(tdfEnv) %*% wt
mah <- StatMatch::mahalanobis.dist(tdfEnv, vc = cov)
diag(mah) <- NA
singular <- which(msing == TRUE | msing0 == TRUE)
outs <- which( sapply(coefs, function(x) { any(x$b > 0 ) } ) == TRUE)
# outs2 <- which( sapply(coefs, function(x) { mean(x$n0) > 5 } ) == TRUE)
undef <- sapply(coefs, function(x) nrow(x) < M )
# remove singular  models from matching  
mah[, singular] <- NA
mah[, undef] <- NA
mah[, outs] <- NA
# mah[, outs2] <- NA

match <- apply(mah, 1, which.min)

matchpred <- list(pred = list(), obs = list())

for(i in 1:N){
print(i)
  obs <- datlist[[i]]
  j <- match[i] # creates an index for a matching site to i
  # glmer out-of-sample Mahalanobis nearest match
  pmatch <- pxl.dist(tpxlcov[[i]], tpxlcov[[j]], coefs, M) 
  n <- nrow(tpxlcov[[i]])
  
  a <- coefs[[j]]$a
  b <- coefs[[j]]$b
  n0 <- coefs[[j]]$n0
  t <- 1:tvec[i]-1
  pars <- matrix(NA, length(pmatch), 3)
  pars[, 1] <- a[pmatch]
  pars[, 2] <- b[pmatch]
  pars[, 3] <- n0[pmatch]
  
  pred.m <- gomp(pars[,1], pars[,2], pars[,3], t)
    
  obs.m <- obs %>% as.data.frame() %>% 
    mutate(cl = tpxlcov[[i]][, paste0("cluster", M)]) %>% 
    group_by(cl) %>% 
    summarize_all(mean) %>% 
    dplyr::select(-cl) %>% as.matrix()
  
  matchpred$pred[[i]] <- pred.m 
  matchpred$obs[[i]] <- obs.m
}
lout[[M]] <- matchpred
# saveRDS(matchpred, file = paste0("nlcd/model_pred/match_pred_", M, ".rds"))

# === matching at the site level only
##################################
#matchsitepred <- list(pred = list(), obs = list())
#
#for(i in 1:N){
#  obs <- datlist[[i]]
#  # out-of-sample Mahalanobis nearest match
#  j <- match[i]
#  a <- mean(coefs[[j]]$a)
#  b <- mean(coefs[[j]]$b)
#  n0 <- mean(coefs[[j]]$n0)
#  t <- 1:tvec[i]-1
#  
#  matpred.m <- matrix(gomp(a, b, n0, t), 
#                      nr = nrow(obs), 
#                      nc = length(t), byrow = T)
#  
#  # glmer out-of-sample, random match
#  l <- sample(1:N, 1, replace = T)
#  a <- mean(coefs[[l]]$a)
#  b <- mean(coefs[[l]]$b)
#  n0 <- mean(coefs[[l]]$n0)
#  t <- 1:tvec[i]-1
#  
#  matpred <- matrix(gomp(a, b, n0, t), 
#                      nr = nrow(obs), 
#                      nc = length(t), byrow = T)
#  matchsitepred$pred[[i]] <- matpred
#  matchsitepred$obs[[i]] <- obs
#}
#
#saveRDS(matchsitepred, file = paste0("nlcd/model_pred/matchsite_pred_", M, ".rds"))

# close the loop at the cluster level
}

# === quick vis
maemat <- matrix(NA, 3, N)
for(i in 1:N){
  # ymu <- matchsite$obs[[i]] %>% 
  #   as.data.frame() %>% 
  #   summarize_all(mean) %>% as.matrix()
  maemat[1,i] = mean(abs(matchsite$obs[[i]][9] - matchsite$pred[[i]][9] ))
  maemat[2,i] = mean(abs(lout[[4]]$obs[[i]][,9] - lout[[4]]$pred[[i]][, 9] ))
  maemat[3,i] = mean(abs(lout[[14]]$obs[[i]][,9] - lout[[14]]$pred[[i]][, 9] ))
}
cols = viridis::viridis(8)
#  --- density plots
# pdf("figures/fig3.pdf", width = 4, height = 4)
plot(density(insample), lwd = 2, xlim = c(0, 20), ylim = c(0, .3),
     xlab = "MAE, % cover", col = cols[1], bty = "n", main = "")
polygon(density(insample), col = col_a(cols[1], .25))
lines(density(maemat[1,]), lwd = 2, col = cols[3])
polygon(density(maemat[1, ]), col = col_a(cols[3], .25))
lines(density(maemat[2, ]), lwd = 2, col = cols[5])
polygon(density(maemat[2, ]), col = col_a(cols[5], .25))
lines(density(maemat[3, ]), lwd = 2, col = cols[7])
polygon(density(maemat[3, ]), col = col_a(cols[7], .25))
# dev.off()
          
plot(maemat[1,] ~ kvec, pch = 19, col = rgb(0, 1, 0, .2), 
     ylab = "MAE, % cover", xlab = "Mean pre-disturbance cover, %", 
     ylim = c(0, 25))
points(maecl[1, ] ~ kvec, pch = 19, col = rgb(0, .5, .5, .2))
points(maecl[2, ] ~ kvec, pch = 19, col = rgb(0, 1, 1, .2))
abline(0, 1, lwd = 3)
points(maenull ~ kvec, pch = 19, col = rgb(.1, 1, 1, .2))

bind_cols(maenull, mae0, maecl[1,], maecl[2,]) %>% 
  pivot_longer(cols = 1:ncol(.)) %>% 
  group_by(name) %>% 
  summarize_all(mean) %>% 
  ungroup() %>% 
  ggplot(aes(x = name, y = value)) +
  geom_point() +
  theme_bw()



# === try weighted mahalanobis distance for matching
library("ParBayesianOptimization")

modout <- readRDS("~/Downloads/modout_1.rds")
datlist <- readRDS("~/Downloads/datlist_1.rds")

coefs <- modout$coefs
tvec <- modout$tvec

N <- length(coefs)

# === optimize
bounds <- list(# dem = c(0.01, 1), demsd = c(0.01, 1), 
               # chili = c(0.01, 1), chilisd = c(0.01, 1))
               ppt = c(0.1, 15), tmax = c(0.1, 15), tmin = c(0.1, 15), 
               sm1 = c(0.1, 15), sm2 = c(0.1, 15), sm3 = c(0.1, 15))
initGrid <- data.frame(# dem = c(seq(.01, 1, l = 15)), 
                       # demsd = c(seq(.01, 1, l = 15)),
                       # chili = c(seq(.01, 1, l = 15)),
                       # chilisd = c(seq(.01, 1, l = 15)))
                       ppt = c(seq(1, 10, l = 15)),
                       tmax = c(seq(1, 10, l = 15)),
                       tmin = c(seq(1, 10, l = 15)),
                       sm1 = c(seq(1, 10, l = 15)),
                       sm2 = c(seq(1, 10, l = 15)),
                       sm3 = c(seq(1, 10, l = 15)))
fbo <- function(#dem, demsd, chili, chilisd)
                   ppt, tmax, tmin,
                   sm1, sm2, sm3)
  {

  # wt <- diag( c(dem, demsd, chili, chilisd, ppt, tmax, tmin, sm1, sm2, sm3) )
  # wt <- diag( c(dem, demsd, chili, chilisd, rep(1, 6) ) )
  wt <- diag( c( rep(1, 4), ppt, tmax, tmin, sm1, sm2, sm3) )
  cov <- var((tdfEnv)) %*% wt
  mah <- StatMatch::mahalanobis.dist((tdfEnv), vc = cov)
  diag(mah) <- NA
  outs <- which( sapply(coefs, function(x) { any(x$b > 0) } ) == TRUE)

  mah[, outs] <- NA

  match <- apply(mah, 1, which.min)
  
  matchsite <- list(pred = list(), obs = list())
  for(i in 1:N){
    # print(i)
    obs <- apply(datlist[[i]], 2, mean)
    # out-of-sample Mahalanobis nearest match
    j <- match[i]
    a <- coefs[[j]]$a
    b <- coefs[[j]]$b
    n0 <- coefs[[j]]$n0
    t <- 1:tvec[i]-1
    matpred.m <- matrix(gomp(a, b, n0, t), 
                        nr = 1, 
                        nc = length(t), byrow = T)
    
    matchsite$pred[[i]] <- matpred.m
    matchsite$obs[[i]] <- obs
  }
  
  m <- rep(NA, N)
  for(i in 1:N){
    yhat <- matchsite$pred[[i]]
    m[i] = mean(abs( matchsite$obs[[i]][1:length(yhat)] - yhat ))
  }

  
  err = sd(m) # norm(m, type = "2")
  return(list(Score = -err))
}

opt <- bayesOpt(FUN = fbo, 
                   bounds = bounds, initGrid = initGrid, iters.n = 9,
                   acq = "ei")

pkgs <- c("dplyr", "tidyr")
sapply(pkgs, require, character.only = T)

years <- c(1985:2018)

source("nlcd/helper_fns.R")
tpxlcov <- readRDS("nlcd/data/tpxlcov.rds")
tdfEnv <- readRDS("nlcd/data/tdfEnv_covars.rds")

# --- site level
modout <- readRDS("nlcd/model_output_test/modout_1.rds")
datlist <- readRDS("nlcd/model_output_test/datlist_1.rds")
  
coefs <- modout$coefs
tvec <- modout$tvec
  
N <- length(coefs)

tdfEnv$ninit <- sapply(coefs, function(x) { x$n0 } )
mah <- StatMatch::mahalanobis.dist(tdfEnv)
diag(mah) <- NA
match <- apply(mah, 1, which.min)

matchsite <- list(pred = list(), obs = list())

for(i in 1:N){
  print(i)
  obs <- datlist[[i]]
  # out-of-sample Mahalanobis nearest match
  j <- match[i]
  a <- coefs[[j]]$a
  b <- coefs[[j]]$b
  n0 <- coefs[[j]]$n0
  t <- 1:tvec[i]-1
  matpred.m <- matrix(gomp(a, b, n0, t), 
                      nr = nrow(obs), 
                      nc = length(t), byrow = T)
  
  matchsite$pred[[i]] <- matpred.m
  matchsite$obs[[i]] <- obs
}

saveRDS(matchsite, file = paste0("nlcd/model_pred/match_pred_1.rds"))


# --- site + cluster levels
clvec <- 2:14
for(M in clvec){
cat("---working on -", M, "-cluster---")  
  inm <- paste0("nlcd/model_output_test/modout_", M, ".rds")
  ind <- paste0("nlcd/model_output_test/datlist_", M, ".rds")
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
tdfEnv$ninit <- sapply(coefs, function(x) { mean(x$n0) } )
mah <- StatMatch::mahalanobis.dist(tdfEnv)
diag(mah) <- NA
singular <- which(msing == T | msing0 == T)
undef <- sapply(coefs, function(x) nrow(x) < M )
# remove singular  models from matching  
mah[, singular] <- NA
mah[, undef] <- NA
match <- apply(mah, 1, which.min)

matchpred <- list(pred = list(), predw = list())

for(i in 1:N){
print(i)
  j <- match[i] # creates an index for a matching site to i
  # glmer out-of-sample Mahalanobis nearest match
  pmatch <- pxl.dist(tpxlcov[[i]], tpxlcov[[j]], coefs, M) 
  n <- nrow(tpxlcov[[i]])
  
  a <- coefs[[j]]$a
  b <- coefs[[j]]$b
  n0 <- coefs[[j]]$n0
  t <- 1:tvec[i]-1
  pars <- matrix(NA, n, 3)
  pars[, 1] <- a[pmatch]
  pars[, 2] <- b[pmatch]
  pars[, 3] <- n0[pmatch]
  
  pred.m <- gomp(pars[,1], pars[,2], pars[,3], t)
  
  matchpred$pred[[i]] <- pred.m 
}

saveRDS(matchpred, file = paste0("nlcd/model_pred/match_pred_", M, ".rds"))

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

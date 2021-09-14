pkgs <- c("dplyr", "tidyr")
sapply(pkgs, require, character.only = T)

years <- c(1985:2018)

source("helper_fns.R")
tpxlcov <- readRDS("data/tpxlcov.rds")
tdfEnv <- readRDS("data/tdfEnv_covars.rds")

clvec <- 3:15
for(m in clvec){
  M <- m
  
  inm <- paste0("nlcd/model_output/modout_", M, ".rds")
  ind <- paste0("nlcd/model_output/datlist_", M, ".rds")
  modout <- readRDS(inm)
  datlist <- readRDS(ind)
  # modout <- readRDS("~/Downloads/modout_4.rds")
  # datlist <- readRDS("~/Downloads/datlist_4.rds")
  
  coefs <- modout$coefs
  tvec <- modout$tvec
  
  N <- length(coefs)
  
# === In-Sample predictions 
##################################
insample <- list(pred = list(), 
                 obs = list())

for(j in 1:N) {
  a <- coefs[[j]]$a
  b <- coefs[[j]]$b
  n0 <- coefs[[j]]$n0
  t <- 1:tvec[j]-1
  matpred <- gomp(a, b, n0, t) %>% 
    as.data.frame() %>% 
    mutate(cl = row_number()) 
  
  jmat <- datlist[[j]] %>%
    as.data.frame() %>%
    mutate(cluster = tpxlcov[[j]][, paste0("cluster", M)]) %>% 
    left_join(matpred, by = c("cluster" = "cl"))
  
  obs <- jmat[, 1:length(t)] %>% as.matrix()
  pred <- jmat %>% select(starts_with("V")) %>% as.matrix()
  insample[[1]][[j]] <- pred
  insample[[2]][[j]] <- obs
}

saveRDS(insample, file = paste0("insample_pred_", M, ".rds"))
##################################


# === Matching at site and pixel levels (pixel-to-cluster)
##################################
mah <- StatMatch::mahalanobis.dist(tdfEnv[,])
diag(mah) <- NA
match <- apply(mah, 1, which.min)
mahdist <- apply(mah, 1, min, na.rm = T)

matchpred <- list(pred = list(), 
                 obs = list(),
                 pred_rand = list())

for(i in 1:N){
  # glmer out-of-sample Mahalanobis nearest match
  j <- match[i] # creates an index for a matching site to i
  a <- coefs[[j]]$a
  b <- coefs[[j]]$b
  n0 <- coefs[[j]]$n0
  t <- 1:tvec[i]-1
  
  matpred <- gomp(a, b, n0, t) %>% 
    as.data.frame() %>%
    mutate(cl = row_number())
  
  pred.m <- pxl.dist(tpxlcov[[i]], tpxlcov[[j]], M) %>% 
    select(cluster) %>% 
    left_join(matpred, by = c("cluster" = "cl")) %>%
    select(starts_with("V"))
  
  # glmer out-of-sample, random match
  l <- sample(1:N, 1, replace = T)
  a <- coefs[[l]]$a
  b <- coefs[[l]]$b
  n0 <- coefs[[l]]$n0
  t <- 1:tvec[i]-1
  
  matpred.r <- gomp(a, b, n0, t) %>% 
    as.data.frame() %>% 
    slice(sample(1:n())) %>%
    mutate(cl = row_number()) 
  
  pred.r <- tpxlcov[[i]] %>%
    mutate(cluster = sample(1:M, n(), replace = T)) %>% 
    select(cluster) %>% 
    left_join(matpred.r, by = c("cluster" = "cl")) %>% 
    select(-cluster)
  
  obs <- datlist[[i]] 
  
  matchpred$pred[[j]] <- pred.m %>% as.matrix()
  matchpred$obs[[j]] <- obs
  matchpred$pred_rand[[j]] <- pred.r %>% as.matrix()
}

saveRDS(matchpred, file = paste0("match_pred_", M, ".rds"))

# === matching at the site level only
##################################
matchsitepred <- list(pred = list(), 
                       obs = list(),
                       pred_rand = list())
for(i in 1:N){
  obs <- datlist[[i]]
  # out-of-sample Mahalanobis nearest match
  j <- match[i]
  a <- mean(coefs[[j]]$a)
  b <- mean(coefs[[j]]$b)
  n0 <- mean(coefs[[j]]$n0)
  t <- 1:tvec[i]-1
  
  matpred.m <- matrix(gomp(a, b, n0, t), 
                      nr = nrow(obs), 
                      nc = length(t), byrow = T)
  
  # glmer out-of-sample, random match
  l <- sample(1:N, 1, replace = T)
  a <- mean(coefs[[l]]$a)
  b <- mean(coefs[[l]]$b)
  n0 <- mean(coefs[[l]]$n0)
  t <- 1:tvec[i]-1
  
  matpred.r <- matrix(gomp(a, b, n0, t), 
                      nr = nrow(obs), 
                      nc = length(t), byrow = T)
  
  matchsitepred$pred[[j]] <- matpred.m
  matchsitepred$obs[[j]] <- obs
  matchsitepred$pred_rand[[j]] <- matpred.r
}

saveRDS(matchsitepred, file = paste0("matchsite_pred_", M, ".rds"))

# close the loop at the cluster level
}
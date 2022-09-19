pkgs <- c("dplyr", "tidyr")
sapply(pkgs, require, character.only = T)

years <- c(1985:2018)

source("helper_fns.R")
# === load data
tpxlcov <- readRDS("data/tpxlcov.rds")
tdfEnv <- readRDS("data/tdfEnv_covars.rds")
tfires <- readRDS("data/tfires.rds")
tsage <- readRDS("data/tsage.rds")
txys <- readRDS("data/txys.rds")

k <- sapply(tpxlcov, function(x) { mean(x$prefire) } )

# --- wildfire level
modout <- readRDS("outputs/models/smodout_1.rds")
datlist <- readRDS("outputs/models/datlist_1.rds")

coefs <- modout$coefs
tvec <- modout$tvec
  
N <- length(coefs)

tdfEnv <- select(tdfEnv, tmax, tmin, dem_mean, chili_mean, sm03mean, sm03sd, ppt)
mah <- mah.dist(tdfEnv)
diag(mah) <- NA
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
saveRDS(matchsite, file = paste0("outputs/matchpred.rds"))

# --- wildfire + cluster levels
lout <- list()
clvec <- 3:14
for(M in clvec){
cat("---working on -", M, "-cluster---")  
  inm <- paste0("outputs/models/modout_", M, ".rds")
  ind <- paste0("outputs/models/datlist_", M, ".rds")
  
  modout <- readRDS(inm)
  datlist <- readRDS(ind)
  
  coefs <- modout$coefs
  tvec <- modout$tvec
  
  msing <- unlist(modout$msing)
  msing0 <- unlist(modout$msing0)
  
  N <- length(coefs)

# === Matching at site and pixel levels (pixel-to-cluster)
mah <- StatMatch::mahalanobis.dist(tdfEnv)
diag(mah) <- NA
singular <- which(msing == TRUE | msing0 == TRUE)
outs <- which( sapply(coefs, function(x) { any(x$b > 0 ) } ) == TRUE)
undef <- sapply(coefs, function(x) nrow(x) < M )
# remove singular/...  models from matching  
mah[, singular] <- NA
mah[, undef] <- NA
mah[, outs] <- NA

match <- apply(mah, 1, which.min)

matchpred <- list(pred = list(), obs = list())

for(i in 1:N){
print(i)
  obs <- datlist[[i]]
  j <- match[i] # creates an index for a matching site to i
  # glmer out-of-sample Mahalanobis nearest match
  if(i == j){
    pmatch <- 1:M
  } else {
    pmatch <- pxl.dist(tpxlcov[[i]], tpxlcov[[j]], coefs, M)   
  }
  

  n <- nrow(tpxlcov[[i]])
  
  a <- coefs[[j]]$a
  b <- coefs[[j]]$b
  n0 <- coefs[[j]]$n0
  t <- 1:tvec[i]-1
  pars <- matrix(NA, M, 3)
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
}

saveRDS(lout, file = "outputs/matchsitecluster.rds")



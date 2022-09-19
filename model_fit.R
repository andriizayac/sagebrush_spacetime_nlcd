pkgs <- c("matrixcalc", "dplyr", "tidyr", "lme4", "raster", "spdplyr")
sapply(pkgs, require, character.only = T)

years <- c(1985:2018)

source("helper_fns.R")
# ====================================
# see pxlmatching.R 
tfires <- readRDS("data/tfires.rds")
tsage <- readRDS("data/tsage.rds")
tpxlcov <- readRDS("data/tpxlcov.rds")
tdfEnv <- readRDS("data/tdfEnv_covars.rds")

kvec <- sapply(tpxlcov, function(x){ mean(x[['prefire']]) })

N <- length(kvec)

# === create a directory to store output
dir.create("outputs/models/", showWarnings = FALSE)

# === glm: wildfire level

mlist <- list()
mlist0 <- list()
y <- list()
yhat <- list()
datlist <- list()
for(i in 1:N){
  # estimate growth and K parameters
  dat <- glm.dat(tsage, tfires, tpxlcov, i, 2) # cluster number is a filler here as it is not used in the model
  temp <- glm(y ~ 1 + log(x), offset = log(x), family = poisson(link = "log"), data = dat) #
  
  # estimate initial population size based the first 5 years of data post-fire
  dat.init <- glm.dat.init(tsage, tfires, tpxlcov, i, 2)
  temp_n0 <- glm(cover ~ 1 + t, family = poisson(link = "log"), data = dat.init) # 
  
  # store models
  mlist[[i]] = temp
  mlist0[[i]] = temp_n0
  
  y[[i]] = dat
  yhat[[i]] = exp(predict(temp))
  
  datlist[[i]] = dat
}

saveRDS(list(mlist = mlist,
             mlist0 = mlist0,
             y = y, yhat = yhat),
        file = paste0("outputs/models/modout_", 1,".rds"))

saveRDS(datlist, paste0("outputs/models/datlist_", 1, ".rds"))

# === glmer: cluster-level models

for(k in 2:15){
  
  mlist <- list()
  mlist0 <- list()
  y <- list()
  yhat <- list()
  datlist <- list()
  for(i in 1:N){
    
    # estimate growth and K parameters
    dat <- glm.dat(tsage, tfires, tpxlcov, i, k)
    temp <- glmer(y ~ -1 + (1+x|cl), offset = log(x), family = poisson(link = "log"), data = dat) #
    
    # estimate initial population size based the first 5 years of data post-fire
    dat.init <- glm.dat.init(tsage, tfires, tpxlcov, i, k)
    temp_n0 <- glmer(cover ~ -1 + (1+t|cl), family = poisson(link = "log"), data = dat.init) # 
    
    # store models
    mlist[[i]] = temp
    mlist0[[i]] = temp_n0
    
    y[[i]] = dat
    yhat[[i]] = exp(predict(temp))
  }
  
  saveRDS(list(mlist = mlist,
               mlist0 = mlist0,
               y = y, yhat = yhat),
          file = paste0("outputs/models/modout_", k,".rds"))
  saveRDS(datlist, paste0("outputs/models/datlist_", k, ".rds"))
}



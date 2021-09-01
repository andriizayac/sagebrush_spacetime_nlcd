pkgs <- c("matrixcalc", "dplyr", "tidyr", "lme4", "raster", "spdplyr")
sapply(pkgs, require, character.only = T)

years <- c(1985:2018)

tfires <- readRDS("data/tfires.rds")
tsage <- readRDS("data/tsage.rds")
tpxlcov <- readRDS("data/tpxlcov.rds")
tdfEnv <- readRDS("data/tdfEnv_covars.rds")

kvec <- sapply(tpxlcov, function(x){ mean(x[['prefire']]) })
# ====================================
# see pxlmatching.R 
tfires <- readRDS("data/tfires.rds")
tsage <- readRDS("data/tsage.rds")
tpxlcov <- readRDS("data/tpxlcov.rds")
tdfEnv <- readRDS("data/tdfEnv_covars.rds")

kvec <- sapply(tpxlcov, function(x){ mean(x[['prefire']]) })

# --- parametrize k prior
# remotes::install_github("aursiber/fitdistrplus")
library(fitdistrplus)
klist <- list()
for(i in 1:length(tsage)) { klist[[i]] = tpxlcov[[i]]$prefire }
kprior <- unlist(klist)
fgamma <- fitdist( kprior[kprior > 0], "gamma" )

# === lmer models and predictions ####
# === glmer 
mlist.glm <- list()
mlist0.glm <- list()
coefs <- list()
datlist.glm <- list()
tvec <- rep(NA, length(kvec))

for(i in 1:3){ #length(tsage)
  print(i)
  
  # estimate growth and K parameters
  dat <- glm.dat(tsage, tfires, tpxlcov, i, 3)
  temp <- lmer(y ~ (1|cl) + (0+x|cl), REML = TRUE, data = dat) #
  
  # estimate initial population size based the first 5 years of data post-fire
  dat.init <- glm.dat.init(tsage, tfires, tpxlcov, i, 3)
  temp_n0 <- glmer(cover ~ (1|cl) + (0+t|cl), family = Gamma(link = "log"), data = dat.init) # 
  
  # store models
  mlist.glm[[i]] = temp
  mlist0.glm[[i]] = temp_n0
  
  datlist.glm[[i]] <- data.matrix(tsage[[i]][,c(tfires$FireYer[i]-1984):33])
  tvec[i] <- ncol( datlist.glm[[i]])
  coefs[[i]] = data.frame(a = coef(temp)$cl[,2], b = coef(temp)$cl[,1], n0 = exp(coef(temp_n0)$cl[,2]))
}

sapply(coefs, function(x) { sum(x$n0 < 0) } ) # check for negative N0


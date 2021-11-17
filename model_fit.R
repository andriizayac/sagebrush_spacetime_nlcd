pkgs <- c("matrixcalc", "dplyr", "tidyr", "lme4", "raster", "spdplyr")
sapply(pkgs, require, character.only = T)

years <- c(1985:2018)

source("nlcd/helper_fns.R")
# ====================================
# see pxlmatching.R 
tfires <- readRDS("nlcd/data/tfires.rds")
tsage <- readRDS("nlcd/data/tsage.rds")
tpxlcov <- readRDS("nlcd/data/tpxlcov.rds")
tdfEnv <- readRDS("nlcd/data/tdfEnv_covars.rds")

kvec <- sapply(tpxlcov, function(x){ mean(x[['prefire']]) })

N <- length(kvec)

# === lmer models and predictions ####

mlist <- list()
mlist0 <- list()
y <- list()
yhat <- list()
for(i in 1:N){
  	# estimate growth and K parameters
	dat <- glm.dat(tsage, tfires, tpxlcov, i, 2) # cluster number is a filler here as it is not used in the model
  	temp <- glm(y ~ 1 + x, offset = log(x), family = poisson(link = "log"), data = dat) #

  	# estimate initial population size based the first 5 years of data post-fire
  	dat.init <- glm.dat.init(tsage, tfires, tpxlcov, i, 2)
  	temp_n0 <- glm(cover ~ 1 + t, family = poisson(link = "log"), data = dat.init) # 
  
  	# store models
  	mlist[[i]] = temp
  	mlist0[[i]] = temp_n0
  	
	y[[i]] = dat
	yhat[[i]] = exp(predict(temp))
}

saveRDS(list(mlist = mlist, 
		mlist0 = mlist0,
		y = y, yhat = yhat), 
        	file = paste0("nlcd/models/modout_", 1,".rds"))
# === glmer 
#
#for(k in 2:15){
#
#	mlist <- list()
#	mlist0 <- list()
#	y <- list()
#	yhat <- list()
#	for(i in 1:N){
#
#  		# estimate growth and K parameters
#  		dat <- glm.dat(tsage, tfires, tpxlcov, i, k)
#  		temp <- glmer(y ~ (1|cl) + (0+x|cl), offset = log(x), family = poisson(link = "log"), data = dat) #
#  
#  		# estimate initial population size based the first 5 years of data post-fire
#  		dat.init <- glm.dat.init(tsage, tfires, tpxlcov, i, k)
#  		temp_n0 <- glmer(cover ~ (1|cl) + (0+t|cl), family = poisson(link = "log"), data = dat.init) # 
#  
#  		# store models
#  		mlist[[i]] = temp
#  		mlist0[[i]] = temp_n0
#  		
#		y[[i]] = dat
#		yhat[[i]] = exp(predict(temp))
#	}
#
#	saveRDS(list(mlist = mlist, 
#			mlist0 = mlist0, 
#             		y = y, yhat = yhat), 
#        	file = paste0("nlcd/models/modout_", k,".rds"))
#}
#


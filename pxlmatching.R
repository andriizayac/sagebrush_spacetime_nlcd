# === prior: raster_processing.R

# === libraries + paths
pkgs <- c("cluster", "raster", "rgdal", "dplyr", "sf")
sapply(pkgs, require, character.only = T)

years <- c(1985:2018)

# === import data: these data represent wildfires that meet 
# the following criteria: single fire, no treatment
fires <- readRDS("data/fpolygons.rds")
fires$FireYer <- as.numeric(fires$FireYer)
sage <- readRDS("data/sagelist.rds")
dfEnv <- readRDS("data/dfEnv_covars.rds")
pxldfEnv <- readRDS("data/pxlcovlist.rds")

N <- nrow(fires)

for(i in 1:3){
  for(j in 1:N) {
    r <- projectRaster(pxldfEnv[[i]][[j]], sage[[j]][[1]])
    sage[[j]] <- addLayer(sage[[j]], r)
  }
}


sagedf <- list()
covdf <- list()
xys <- list()
for(i in 1:N){
  dat <- as.data.frame(sage[[i]], xy = TRUE, na.rm = T)
  xys[[i]] <- dat[, 1:2]
  sagedf[[i]] <- dat[, 3:35]
  covdf[[i]] <- dat[, 36:38]
}


# === out of  the larger pool of sites select those that
# after the fire have a decrease in average cover > 1%

ssize <- sapply(sagedf, nrow) # range of sample sizes

maxdiff <- rep(NA, N)
for(i in 1:N){
  l <- as.numeric(apply(sagedf[[i]], 2, function(x) { mean(x[1:23]) } ))
  maxdiff[i] <- min(diff(l))
}

diffins <- which(sapply(maxdiff, function(x) { -x > 1 }))

subid <- diffins


# === add pre-disturbance covariates to cov df
# assigns average stability to the pixels with no variation.

tfires <- fires[subid, ]
tsage <- sagedf[subid]
tdfEnv <- dfEnv[subid, ]
tpxlcov <- covdf[subid]
txys <- xys[subid]

# remove NORTH BLUE MOUNTAIN fire as there is no variation 
outid <- which(tfires$FireNam == "NORTH BLUE MOUNTAIN")
tfires <- tfires[-outid, ]
tsage <- tsage[-outid]
tdfEnv <- tdfEnv[-outid, ]
tpxlcov <- tpxlcov[-outid]
txys <- txys[-outid]

for(i in 1:length(tpxlcov)){
  var1 <- apply(tsage[[i]][, 1:(tfires$FireYer[i]-1985)], 1, mean)
  var2 <- as.numeric(apply(tsage[[i]][, 1:(tfires$FireYer[i]-1985)], 1, function(x) { mean(x)/sd(x) } ))
  tpxlcov[[i]]$prefire <- var1
  stabmean  <- mean(var2[is.finite(var2)], na.rm = T)
  tpxlcov[[i]]$stab <- ifelse(!is.finite(var2), stabmean, var2)
}

# === import/calculate pixel-level covariates

# === apply kmeans clustering
library(cluster)
library(factoextra)

# add cluster 
for(M in 2:15) {
for(i in 1:length(tpxlcov)){
  set.seed(123)
  # https://stackoverflow.com/questions/16274788/k-means-and-mahalanobis-distance 
  # cov(X) = R'R
  # y = XR^-1
  X <- as.matrix(tpxlcov[[i]][,1:5])
  # Re-scale the data
  # C <- chol(var(X))
  # y <- X %*% qr.solve(C)
  y = scale(X)
  k2 <- kmeans(y, centers = M, iter.max = 200, algorithm = "MacQueen")
  tpxlcov[[i]][, paste0("cluster", M)] <- as.numeric(k2$cluster)
}
}

# export data: these are the data that were used in the analysis
saveRDS(tfires, "data/tfires.rds")
saveRDS(tpxlcov, "data/tpxlcov.rds")
saveRDS(tsage, "data/tsage.rds")
saveRDS(tdfEnv, "data/tdfEnv_covars.rds")
saveRDS(txys, "data/txys.rds")
# === next: model_fit.R

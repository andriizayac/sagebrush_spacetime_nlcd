# === libraries + paths
pkgs <- c("cluster", "raster", "rgdal", "spdplyr")
sapply(pkgs, require, character.only = T)

years <- c(1985:2018)
path <- "~/Desktop/SSM_PDE/nlcd_subset/"

# === import data
fires <- readRDS(paste0(path, "fpolygons.rds"))
fires$FireYer <- as.numeric(fires$FireYer)
sage <- readRDS(paste0(path,"sagelist.rds"))
dfEnv <- readRDS(paste0(path,"dfEnv_covars.rds"))
pxldfEnv <- readRDS(paste0(path,"pxlcovlist.rds"))

N <- nrow(fires)

for(i in 1:3){
  for(j in 1:N) {
    r <- projectRaster(pxldfEnv[[i]][[j]], sage[[j]][[1]])
    sage[[j]] <- addLayer(sage[[j]], r)
  }
}


sagedf <- list()
covdf <- list()
for(i in 1:N){
  dat <- as.data.frame(sage[[i]], xy = FALSE, na.rm = T)
  sagedf[[i]] <- dat[, 1:33]
  covdf[[i]] <- dat[, 34:36]
}


# === subset 12 test fires
ssize <- sapply(sagedf, nrow) # range of sample sizes

maxdiff <- rep(NA, N)
for(i in 1:N){
  l <- as.numeric(apply(sagedf[[i]], 2, function(x){mean(x[1:24])}))
  maxdiff[i] <- min(diff(l))
}

diffins <- which(maxdiff < quantile(maxdiff, .33, na.rm = T))
sizeins <- which(ssize > 200)


subid <- sizeins[sizeins %in% diffins]
# test set, manually picked

par(mfrow = c(4,3), mar = c(1,2,1,1))
for(i in 1:12){
  j = subid[i]
  matplot(t(sagedf[[j]]), type = "l", col = rgb(.5, 0, .5, .05), main = j)
  abline(v = fires$FireYer[j]-1984, lwd = 2)
}

# === add pre-disturbance covariates to cov df
tfires <- fires[subid, ]
tsage <- sagedf[subid]
tdfEnv <- dfEnv[subid, ]
tpxlcov <- covdf[subid]

for(i in 1:length(subid)){
  var1 <- apply(tsage[[i]][,1:(tfires$FireYer[i]-1985)], 1, mean)
  var2 <- as.numeric(apply(tsage[[i]][,1:(tfires$FireYer[i]-1985)], 1, function(x){ mean(x)/sd(x) }))
  tpxlcov[[i]]$prefire <- var1
  stabmean  <- mean(var2[is.finite(var2)], na.rm = T)
  tpxlcov[[i]]$stab <- ifelse(!is.finite(var2), stabmean, var2)
}
# assigns average stability to the pixles with no variation.

# === import/calculate pixel-level covariates

# === apply kmeans clustering
library(cluster)
library(factoextra)

# add cluster 
for(i in 1:length(subid)){
  set.seed(123)
  # https://stackoverflow.com/questions/16274788/k-means-and-mahalanobis-distance 
  # cov(X) = R'R
  # y = XR^-1
  X <- as.matrix(tpxlcov[[i]])
  # Rescale the data
  C <- chol(var(X))
  y <- X %*% solve(C)
  k2 <- kmeans(y, centers = 10)
  tpxlcov[[i]]$cluster <- as.numeric(k2$cluster)
}

# export data
saveRDS(tfires, paste0(path, "tfires.rds"))
saveRDS(tpxlcov, paste0(path, "tpxlcov.rds"))
saveRDS(tsage, paste0(path, "tsage.rds"))
saveRDS(tdfEnv, paste0(path, "tdfEnv_covars.rds"))







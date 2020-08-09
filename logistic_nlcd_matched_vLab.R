library(rstan)
library(raster)
library(rgdal)
library(leaflet)
library(viridis)
library(reshape2)

# === allow multi-core stan sa,mpling
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


years <- c(1985:2017)
path <- "C:/Users/CaughlinLab/Desktop/Landsat_eros/"

sage <- readRDS(paste0(path, "usgs_sagebrush/sagelist_finite.rds"))
dfspat <- readRDS(paste0(path, "mtch_polygons.rds"))

df <- spTransform(dfspat, sage[[1]][[1]]@crs)
plot(df[c(2, 43), ])
plot(sage[[df$idvar[2]]][[1]], add = T)
plot(sage[[df$idvar[43]]][[1]], add = T)

#dfsubset = df[df$ha_clipped < 50 & df$Fire_Year > 1990,]
#sage1 = subset(sage, df$idvar %in% dfsubset$idvar)

plot(sage[[1]][[1]])
plot(dfspat[1,],add=T)
################## find the fire year based on the trajectories; Create offset variable
n <- nrow(dfspat)
offset <- rep(NA, n)
yrs <- dfspat$Fire_Year - 1984
sagemeans <- matrix(NA, nrow = n, ncol = 33)
for(i in 1:n){
  sagemeans[i, ] <- cellStats(sage[[i]], mean)
  f <- which.min(diff(sagemeans[i, ])) + 1
 if (yrs[i] == f - 1) {
    offset[i] = 1
  }  else if (yrs[i] == f - 2) {
    offset[i] = 2
  } else if (yrs[i] == f + 1) {
    offset[i] = -1
  } else {
    offset[i] = 0
  }
}
offset[250] <- -1
offset[25] <- 1
ind <- which(offset == 0)
ind <- 1:nrow(dfspat)

k=1
fig()
####################################################################

############################################## dopar
library(foreach)
library(doParallel)
cl=makeCluster(20)
registerDoParallel(cl)
#modlist = 
foreach(i = 1:100, .packages=c('raster', 'rstan')) %dopar% {
  # alist <- dataprep(n = i, pol = dfspat, rast = sage, year = years, yroffset = offset, trsubset = .25,
  #                predsubset = 1, threshold = TRUE, thresholdval = .1)
  # saveRDS(alist, paste0(path, "models_stan_eros/datasets_models_eros/dat", i, ".rds"))
  alist <- readRDS(paste0(path,"models_stan_eros/datasets_models_eros/dat",i,".rds"))
  temp <- stan(file = paste0(path, "rdiff_discrete_landsat_non_spat_statespace.stan"),
              data = alist, iter = 500, warmup = 450, chains = 3, seed = 125,
              pars = c("Z", "r", "r_mu", "r_sigma", "k", "k_mu", "k_phi", "k_sigma", "eta",
                       "y_pred", "y_pred01", "z_pred", "z_init", "gamma", "alpha", "alpha_pred", "sigma_alpha_pred"),
              save_warmup = FALSE)
  saveRDS(temp, paste0(path, "models_stan_eros/mod", i, ".rds"))
}
stopCluster(cl)
#####################################################
# === explore the fit models
modlist <- list()
datlist <- list()
k_mean <- rep(NA, 25)
diverge <- rep(NA, 25)
for(i in 1:25) {
  print(i)
  temp <- readRDS(paste0(path, "/models_stan_eros/", "mod", i, ".rds"))
  diverge[i] <- get_num_divergent(temp)
  modlist[[i]] <- rstan::extract(temp)
  datlist[[i]] <- readRDS(paste0(path, "/models_stan_eros/datasets_models_eros/", "dat", i, ".rds"))
  k_mean[i] <- mean(datlist[[i]]$y_k_pred)
}


# --- mean absolute error
plotTraj(modlist[[1]], datlist[[1]])
#matplot(t(as.matrix(alist$dftr)),type="l",add=F,col=rgb(0,0,0,.25))
plot(t(yp[, ]) ~ jitter(as.matrix(datlist$dftr)), pch = 3, xlab = "NLCD", ylab = "Predicted", col = rgb(0,0,0,.5))
abline(0,1,col="red",lwd=3)

y_pred <- with(post, y_pred)
err <- rep(NA, dim(y_pred)[1])

for(i in 1:length(err)){err[i] = mae(t(y_pred[i,,]),as.matrix(datlist$dftest))}
plot(density(err))
polygon(density(err), col = "blue")

### ---set up error data frame
postsample <- length(with(modlist[[1]], alpha))
err <- matrix(NA, postsample, 25)
col <- viridis(25)
for(j in 1:ncol(err)){
  print(j)
  post <- modlist[[j]]
  y_pred <- with(post, y_pred)
  for(i in 1:dim(y_pred)[1]){
    err[i, j] <- mae(y_pred[i, , ], as.matrix(datlist[[j]]$dftest))
  }
}

plot(density(err[, 1]), col = col[1], xlim=c(0, 10), ylim = c(0, 15), main = "Out-of-sample error", xlab="MAE")
for(i in 1:ncol(err)){polygon(density(err[, i]), col = col[i])}
#legend("topright", legend = c("simple logistic", "+varying k, +proc err", "spatial"), fill=c(col[c(1,6)],"red"))

# === explore variation in predictive capacity
mahvar <- rep(NA, 25)
for(i in 1:length(mahvar)) { 
  mahvar[i] = mahdist[ids[i], idsRef[i]]
}

cbind(dfEnv[1:25, ], k_mean, mahvar, diverge, t(err)) %>% melt(id.vars = c(names(dfEnv), "k_mean", "mahvar", "diverge")) %>% 
  #ggplot(aes(y = value, x = as.factor(idvar))) + geom_boxplot() +
  ggplot(aes(y = value, x = diverge)) + geom_point() +
  labs(y = "MAE")



# ==== plot fire polygons with leaflet ####
poly <- spTransform(df, CRS("+proj=longlat +datum=WGS84"))
ids <- 1:25
idsRef <- poly$closeMatch[1:25]

r1 <- sage[[ids[1]]][[1]]
r2 <- sage[[ids[1]]][[1]]
#crs(r1) = sp::CRS("+init=epsg:3857")

pal <- colorNumeric(viridis(5), values(r1),
                    na.color = "transparent")

leaflet() %>% 
  addPolygons(data = poly[ids, ], color = "blue", fillOpacity = .9) %>%
  addPolygons(data = poly[idsRef, ], color = "red", fillOpacity = .9) %>%
  addTiles(group = "OSM",
                       options = providerTileOptions(minZoom = .1, maxZoom = 100)) %>%
  addProviderTiles("Esri.WorldTopoMap",    
                   group = "Topo") 
  addRasterImage(r1, colors = pal, opacity = 0.8,project = TRUE) %>%
  addRasterImage(r2, colors = pal, opacity = 0.8,project = TRUE) %>%
  addLegend(pal = pal, values = values(r1),
            title = "Sagebrush cover, %")

#########################


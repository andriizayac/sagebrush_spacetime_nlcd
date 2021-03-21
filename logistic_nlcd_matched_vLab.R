library(brms)
library(raster)
library(rgdal)
library(leaflet)
library(viridis)
library(reshape2)
library(dplyr)
library(ggplot2)

# === allow multi-core stan sa,mpling
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

years <- c(1985:2017)
path <- "C:/Users/CaughlinLab/Desktop/Landsat_eros/"

# prepare intput data for non-spatial stan model

years <- c(1985:2018)
path <- "~/Desktop/SSM_PDE/nlcd_subset/"
fires <- readOGR(paste0(path,"GB_wildfire_subset_1951_2018_ct_1_Rx0_blm_epsg4326.shp"))[1:10,]

sage <- readRDS(paste0(path,"sagelist_subs.rds"))
sagecrs <- sage[[1]]@crs

df <- spTransform(fires, sagecrs)
plot(df[c(4),])
plot(sage[[df$idvar[1]]][[1]],add=T)

#dfsubset = df[df$ha_clipped < 50 & df$Fire_Year > 1990,]
#sage1 = subset(sage, df$idvar %in% dfsubset$idvar)

mat <- as.data.frame(sage[[1]], xy = FALSE, na.rm = TRUE)
matplot(t(mat), type = "l", col = viridis(10), ylab = "Cover, %", xlab = "Time, years")
abline(v = 22)
# === find the fire year based on the trajectories; Create offset variable
# n <- nrow(dfspat)
# offset <- rep(NA, n)
# yrs <- dfspat$Fire_Year - 1984
# sagemeans <- matrix(NA, nrow = n, ncol = 33)
# for(i in 1:n){
#   sagemeans[i, ] <- cellStats(sage[[i]], mean)
#   f <- which.min(diff(sagemeans[i, ])) + 1
#  if (yrs[i] == f - 1) {
#     offset[i] = 1
#   }  else if (yrs[i] == f - 2) {
#     offset[i] = 2
#   } else if (yrs[i] == f + 1) {
#     offset[i] = -1
#   } else {
#     offset[i] = 0
#   }
# }
# offset[250] <- -1
# offset[25] <- 1
# ind <- which(offset == 0)
# ind <- 1:nrow(dfspat)
# 
# k=1
# fig()
# =============================================================
# see pxlmatching.R 
tfires <- readRDS(paste0("C:/Users/CaughlinLab/Desktop/Landsat_eros/tfires.rds"))
tsage <- readRDS(paste0("C:/Users/CaughlinLab/Desktop/Landsat_eros/tsage.rds"))
tpxlcov <- readRDS(paste0("C:/Users/CaughlinLab/Desktop/Landsat_eros/tpxlcov.rds"))

############################################## dopar 1
library(foreach)
library(doParallel)
cl=makeCluster(20)
registerDoParallel(cl)

foreach(i = 1:length(tsage), .packages=c('reshape2','brms','matrixcalc')) %dopar% {
  # --- Gompertz Poisson
  mat <- data.matrix(tsage[[i]][,c(tfires$FireYear[i]-1984):33])
  xpr <- vec(mat[,-ncol(mat)])
  dat <- data.frame(y = vec(mat[,-1]), 
                    x = log(ifelse(xpr == 0, .2, xpr)),
                    cl = as.factor(rep(tpxlcov[[i]]$cluster, ncol(mat)-1)))
  
  temp <- brm(y ~ 1 + x + (1|cl) + offset(x), family="poisson", data=dat, chains = 4, iter = 1000, warmup = 800)
  
  saveRDS(temp, paste0(path, "models_stan_eros/gomp_log_", i, ".rds"))
  rm(temp)
}
stopCluster(cl)
#####################################################

############################################## dopar 2
library(foreach)
library(doParallel)
cl=makeCluster(20)
registerDoParallel(cl)

foreach(i = 1:length(tsage), .packages=c('reshape2','brms','matrixcalc')) %dopar% {
  # --- exact Gompertz solution
  # u(t) = k*exp(const*exp(-alpha*t)), where const = log(N0/k)
  mat <- tsage[[i]][,c(tfires$FireYear[i]-1984):33]
  colnames(mat) <- 1:ncol(mat)
  mat$n0 = as.numeric(ifelse(mat[,1] == 0, .2, mat[,1]))
  mat$cl = as.factor(tpxlcov[[i]]$cluster)
  dat <- melt(mat, id.vars = c("n0", "cl"))
  names(dat) <- c("n0", "cl", "t", "cover")
  dat$t <- as.numeric(dat$t)
  dat$cover = ifelse(dat$cover==0, .2, dat$cover)
  
  gompnl <- bf(cover ~ k*exp(log(n0/k)*exp(-alpha*t)),
               alpha ~ 1 + (1|cl), k ~ 1, nl = TRUE)
  nlprior <- c(prior(normal(10, 5), nlpar = "k"),
               prior(normal(0, 1), nlpar = "alpha"))
  temp <- brm(formula = gompnl, data = dat, family = gaussian(),
              prior = nlprior, 
              control = list(adapt_delta = 0.9),
              chains = 4, iter = 1000, warmup = 800)
  
  saveRDS(temp, paste0(path, "models_stan_eros/gomp_", i, ".rds"))
  rm(temp)
}
stopCluster(cl)
#####################################################

############################################## dopar 3
library(foreach)
library(doParallel)
cl=makeCluster(20)
registerDoParallel(cl)

foreach(i = 1:length(tsage), .packages=c('reshape2','brms','matrixcalc')) %dopar% {
  # --- exact Gompertz solution
  # u(t) = k*exp(const*exp(-alpha*t)), where const = log(N0/k)
  mat <- tsage[[i]][,c(tfires$FireYear[i]-1984):33]
  colnames(mat) <- 1:ncol(mat)
  mat$n0 = as.numeric(ifelse(mat[,1] == 0,  0.01, mat[,1]))
  mat$cl = as.factor(tpxlcov[[i]]$cluster)
  dat <- melt(mat, id.vars = c("n0", "cl"))
  names(dat) <- c("n0", "cl", "t", "cover")
  dat$t <- as.numeric(dat$t)
  dat$cover = ifelse(dat$cover==0, 0.01, dat$cover)
  
  lgnl <- bf(cover ~ k*n0*exp(r*t)/(k+n0*(exp(r*t)-1) ),
             r ~ 1 +(1|cl), k ~ 1, nl = TRUE)
  nlprior <- c(prior(normal(10, 5), nlpar = "k"),
               prior(normal(0, 1), nlpar = "r"))
  fit1 <- brm(formula = lgnl, data = dat, family = gaussian(),
              prior = nlprior, 
              control = list(adapt_delta = 0.9),
              chains = 4, iter = 1000, warmup = 800)
  
  saveRDS(temp, paste0(path, "models_stan_eros/gomp_lg", i, ".rds"))
  rm(temp)
}
stopCluster(cl)
#####################################################


# === explore the fit models
modlist <- list()
datlist <- list()
k_mean <- rep(NA, 100)
diverge <- rep(NA, 100)
for(i in 1:100) {
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
err <- matrix(NA, postsample, 100)
col <- viridis(100)
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
mahdist <- readRDS("mahdist_matrix.rds")
dfEnv <- readRDS("dfEnv_covars.rds")
mahvar <- rep(NA, 100)
for(i in 1:length(mahvar)) { 
  mahvar[i] = mahdist[ids[i], idsRef[i]]
}

a <- cbind(dfEnv[1:100, ], ha_clipped = dfspat$ha_clipped[1:100],  k_mean, mahvar, diverge, t(err)) %>% 
  melt(id.vars = c(names(dfEnv), "k_mean", "mahvar", "diverge", "ha_clipped")) %>% 
  group_by(idvar) %>% mutate(lower = quantile(value, .025), upper = quantile(value, .975), mean = mean(value)) %>% ungroup() #%>% 
  ggplot(aes(y = mean, x = k_mean)) + geom_pointrange(aes(ymin = lower, ymax = upper),size = .1) +
  # scale_colour_viridis("Pre-disturbance abundance") +
  labs(y = "MAE", x = "Pre-disturbance cover, %") +geom_abline(intercept = 0, slope = 1) +
  theme_bw()



# ==== plot fire polygons with leaflet ####
poly <- spTransform(df, CRS("+proj=longlat +datum=WGS84"))
ids <- 1:100
idsRef <- poly$closeMatch[1:100]

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

######################### export figures
library(export)
pathexport <- "C:/Users/CaughlinLab/Downloads/"
graph2eps(file  = paste0(pathexport, "mae.eps"), width = 8, height = 6)

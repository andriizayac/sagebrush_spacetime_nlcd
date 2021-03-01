library(rstan)
library(raster)
library(rgdal)
library(leaflet)
library(viridis)

# === allow multi-core stan sa,mpling
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# prepare intput data for non-spatial stan model

years <- c(1985:2017)
path <- "~/Google Drive/Data/nlcd_fires_matched/"
fires <- readOGR("~/Google Drive/GEE_seminar/Wildfires_1870_2015_Great_Basin_SHAPEFILE/Wildfires_1870_2015_Great_Basin.shp")

sage <- readRDS(paste0(path,"sage_subset/sagelist_finite.rds"))
dfspat <- readRDS(paste0(path,"mtch_polygons.rds"))

df <- spTransform(dfspat,sage[[1]][[1]]@crs)
plot(df[c(1),])
plot(sage[[df$idvar[1]]][[1]],add=T)
plot(sage[[df$idvar[43]]][[1]],add=T)

#dfsubset = df[df$ha_clipped < 50 & df$Fire_Year > 1990,]
#sage1 = subset(sage, df$idvar %in% dfsubset$idvar)
i <- 1
mat <- as.data.frame(sage[[3]], xy = FALSE, na.rm = TRUE)[1:1000, ]
matplot(t(mat), type = "l", col = viridis(10), ylab = "Cover, %", xlab = "Time, years")

plot(sage[[1]][[1]])
plot(dfspat[1,],add=T)
################## find the fire year based on the trajectories; Create offset variable
n <- nrow(dfspat)
offset <- rep(NA, n)
yrs <- dfspat$Fire_Year - 1984
sagemeans <- matrix(NA, nrow = n, ncol = 33)
for(i in 1:n){
  sagemeans[i, ] = cellStats(sage[[i]], mean)
  f = which.min(diff(sagemeans[i, ])) + 1
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
# plot the ts of average shrub values with proposed fire time point
k = 9
fig(index = dfspat$closeMatch)
####################################################################

### input data
alist <- dataprep2(n = 1, pol = dfspat, rast = sage, year = years, yroffset = offset, trsubset = .25, 
                   predsubset = 1, threshold = TRUE, thresholdval = .1)
matplot(t(as.matrix(alist$dftr)), type = "l")
# --- fit the model to the tr fire and predict over the test fire
model2a <- stan_model("~/Desktop/CodingClub/stan_models/rdiff_discrete_landsat_non_spat_statespace.stan")
model2 <- stan_model("~/Desktop/CodingClub/stan_models/rdiff_discrete_landsat_non_spat.stan")
mod2a_1 <- sampling(model2a, data = alist, iter = 50, warmup = 45, chains = 3, seed = 125,
             pars = c("Z", "r", "r_mu", "r_sigma", "k", "k_mu", "k_phi", "k_sigma", "eta", 
                      "y_pred", "y_pred01", "z_pred", "z_init", "gamma", "alpha"),
             save_warmup = FALSE)

# --- mean aboslute error
post <- rstan::extract(mod2a_1)
yp <- apply(with(post, Z), c(2, 3), mean)#-matrix(rep(gamma,930),ncol=930)
matplot(yp[,], type = 'l')
#matplot(t(as.matrix(alist$dftr)),type="l",add=F,col=rgb(0,0,0,.25))
plot(t(yp) ~ jitter(as.matrix(alist$dftest)), pch = 3, xlab = "NLCD", ylab = "Predicted", col = rgb(0,0,0,.5))
abline(0, 1, col = "red", lwd = 3)

y_pred <- with(post, y_pred)
err <- rep(NA, dim(y_pred)[1])

for(i in 1:length(err)){ err[i] = mae(t(y_pred[i, , ]), as.matrix(alist$dftest)) }
plot(density(err))
polygon(density(err), col = "blue")
### ---set up error data frame
err <- matrix(NA, dim(y_pred)[1], 8)
col <- viridis(8)
for(j in 1:ncol(err)){
  post <- rstan::extract(get(paste0("m", j)))
  y_pred <- with(post, y_pred)
  for(i in 1:dim(y_pred)[1]){
    err[i,j] <- mae(y_pred[i, , ], mat[test, ])
  }
}
plot(density(err[, 1]), col = col[1], xlim = c(1.5, 5), ylim = c(0, 15), main = "Out-of-sample error", xlab = "MAE")
for(i in 1:ncol(err)){ polygon(density(err[, i]), col = col[i]) }
legend("topright",legend=c("simple logistic","+varying k, +proc err","spatial"),fill=c(col[c(1,6)],"red"))

############################################## dopar
library(foreach)
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
# modlist = 
foreach(i = 1:25, .packages = c('raster','rstan')) %dopar% {
  alist <- dataprep2(n = i, pol = dfspat, rast = sage, year = years, yroffset = offset, trsubset = .25,
                 predsubset = 1,threshold = TRUE,thresholdval = .1)
  # temp =stan(file="~/Desktop/CodingClub/stan_models/rdiff_discrete_landsat_non_spat_statespace.stan",
  #            data=alist,iter=50,warmup=45,chains=3,seed=125,
  #            pars=c("Z","r","r_mu","r_sigma","k","k_mu","k_phi","k_sigma","eta","y_pred","y_pred01","z_pred","z_init","gamma","alpha"),
  #            save_warmup=FALSE)
  saveRDS(alist,paste0("~/Downloads/dat",i,".rds"))
}
stopCluster(cl)

################################## 
res1 <- readRDS("~/Downloads/res1.rds")
ggplot(res1 , aes(y = value, x = as.factor(idvar))) + # geom_pointrange(aes(ymin = lower, ymax = upper), size = .1) +
  geom_boxplot() + 
  # scale_colour_viridis("Pre-disturbance abundance") +
  labs(y = "MAE", x = "ID") +  # geom_abline(intercept = 0, slope = 1) +
  theme_bw() + theme(axis.text.x = element_blank())

res1a <- res1[which(res1$mean < 5 & res1$k_mean > 10), ]

#####################################################
# --- plot fire polygons with leaflet ####
poly <- spTransform(fires, CRS("+proj=longlat +datum=WGS84"))
ids1 <- poly$idvar[1:100]
ids1Ref <- poly$closeMatch[1:100]

r1=sage[[ids1[1]]][[1]]
r2=sage[[ids1[2]]][[1]]
#crs(r1) = sp::CRS("+init=epsg:3857")

pal <- colorNumeric(viridis(5), values(r1),
                    na.color = "transparent")
leaflet(poly[poly$Fire_Name == "Soda",]) %>%
  addPolygons(stroke = T, fillColor = "red", fillOpacity = 1) %>%
  addTiles(group = "OSM",
           options = providerTileOptions(minZoom = .1, maxZoom = 100)) %>%
  addProviderTiles("Esri.WorldTopoMap",    
                   group = "Topo") %>%
  addScaleBar(position = "bottomright")


leaflet() %>% 
  # addPolygons(data = poly, color = "black") %>%
  addPolygons(data = poly[ids1, ], color = "red") %>%
  addPolygons(data = poly[ids1Ref, ], color = "blue") %>%
  addTiles(group = "OSM",
                       options = providerTileOptions(minZoom = .1, maxZoom = 100)) %>%
  addProviderTiles("Esri.WorldTopoMap",    
                   group = "Topo") %>% 
  addScaleBar(position = "bottomright")

  # addRasterImage(r1, colors = pal, opacity = 0.8,project = TRUE) %>%
  # addRasterImage(r2, colors = pal, opacity = 0.8,project = TRUE) %>%
  # addLegend(pal = pal, values = values(r1),
  #           title = "Sagebrush cover, %")


#### Export figures
# ----
library(export)
pathfig <- "~/Downloads/"
graph2eps(file = paste0(pathfig, "figh2.eps"), height = 6, width = 6)

library(htmlwidgets)
saveWidget(m, file = paste0(pathfig, "figmap.html"))
library(mapview)
mapshot(m, file = paste0(pathfig, "figmap1.pdf"))

dat1 <- readRDS("~/Downloads/dat1.rds")
dat9 <- readRDS("~/Downloads/dat9.rds")
mod1 <- readRDS("~/Downloads/mod1.rds")
mod9 <- readRDS("~/Downloads/mod9.rds")


# --- mean aboslute error
post <- rstan::extract(mod9)
yp <- apply(with(post, y_pred), c(2, 3), mean)#-matrix(rep(gamma,930),ncol=930)
matplot(apply(t(refdat),1,mean), type = 'l', col = viridis(10), 
        ylab = "Cover, %", xlab = "Time, years", main = "Training data")
#matplot(t(as.matrix(alist$dftr)),type="l",add=F,col=rgb(0,0,0,.25))
refdat <- as.data.frame(dat9$dftr, xy = FALSE, na.rm = TRUE)
plot(t(yp) ~ jitter(as.matrix(dat9$dftest)), pch = 3, cex = .25, xlab = "NLCD", ylab = "Predicted", col = rgb(0,0,0,.5))
abline(0, 1, col = viridis(10)[5], lwd = 3)

y_pred <- with(post, y_pred)
err <- rep(NA, dim(y_pred)[1])

for(i in 1:length(err)){ err[i] = mae(t(y_pred[i, , ]), as.matrix(alist$dftest)) }
plot(density(err))
polygon(density(err), col = "blue")


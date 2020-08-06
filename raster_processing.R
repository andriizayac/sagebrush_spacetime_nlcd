library(raster)
library(sf)
library(rgdal)
library(foreach)
library(doParallel)
library(dplyr)
library(spdplyr)
################################ This script collects env data and matches sites using Mahalanobis distance metric

# --- import fire  polygons
# --- fires between 40 and 5000 ha in area burned once
path_fires <- "C:/Users/CaughlinLab/Desktop/Landsat_eros/small_fires_burned_once_40_5000/"
file <- "small_fires_burned_once_40_5000_13scenes.shp"
fire_polygons <- readOGR(paste0(path_fires, file))%>%
  filter(Fire_Year > 1986, ha_clipped > 100)

N <- nrow(fire_polygons)

# --- import rasters to crop --- THIS SECTION OF CODE DOES NOT NEED TO BE REPEATED

# --- 1. Digital Elevation Model: USGS SRTM (30m): GB_SRTM.tif
dem_raster <- raster("C:/Users/CaughlinLab/Desktop/Landsat_eros/GB_SRTM.tif")
demlist <- list()
# --- 2a. Sagebrush NLCD layer (30m)
sage_raster <- stack("C:/Users/CaughlinLab/Desktop/Landsat_eros/usgs_sagebrush/gb_13scenes_sage_mosaic.img")

# --- 2b. Shrub NLCD layer (30m)
shrub_raster <- stack("C:/Users/CaughlinLab/Desktop/Landsat_eros/usgs_shrub/gb_13scenes_shrub_mosaic.img")

# --- 3. CHILI: Continuous Hea-Isolation Load Index (10m): CHILI_30m.tif
#chili_raster=raster("C:/Users/CaughlinLab/Desktop/Landsat_eros/CHILI_30m.tif")
chililist=list()

# --- 4. Climate data (PRISM at 90 m resolution)
# annual precipitation (cumulative): prism_pptAnn
#ppt_raster=raster("C:/Users/CaughlinLab/Desktop/Landsat_eros/prism_pptAnn.tif")
pptlist=list()
# January minimum temperature (C): prism_tmin_jan
#tmin_raster=raster("C:/Users/CaughlinLab/Desktop/Landsat_eros/prism_tmin_jan.tif")
tminlist=list()
# July maximum temperature (C): prism_tmax_jul
#tmax_raster=raster("C:/Users/CaughlinLab/Desktop/Landsat_eros/prism_tmax_jul.tif")
tmaxlist=list()

# --- reeproject fire polygons, clip, and mask the rasters
fpolygons=spTransform(fire_polygons, sage_raster@crs)
a=mask(crop(sage_raster,fpolygons[2,]),fpolygons[2,])
plot(a[[1]])
plot(fpolygons[2,],add=TRUE)

for(i in 1:N){
  print(i)
  #sagebrushlist[[i]] = mask(crop(sage_raster,fpolygons[i,]),fpolygons[i,])
  #demlist[[i]] = mask(crop(dem_raster,fpolygons[i,]),fpolygons[i,])
  #chililist[[i]] = mask(crop(chili_raster,fpolygons[i,]),fpolygons[i,])
  #pptlist[[i]] = mask(crop(ppt_raster,fpolygons[i,]),fpolygons[i,])
  #tminlist[[i]] = mask(crop(tmin_raster,fpolygons[i,]),fpolygons[i,])
  #tmaxlist[[i]] = mask(crop(tmax_raster,fpolygons[i,]),fpolygons[i,])
}

# crop shagerush and shrub rasterstacks ##################################
cl = makeCluster(20)
registerDoParallel(cl)

flist = foreach(i=1:N) %dopar% {
  temp = raster::mask(raster::crop(sage_raster,fpolygons[i,]),fpolygons[i,])
  raster::values(temp)[raster::values(temp)==101]=NA
  temp
}
stopCluster(cl)
######################################################
# ---------------------------------------------------------------- PROCEED FROM HERE
path <- "C:/Users/CaughlinLab/Desktop/Landsat_eros/"
### === Import clipped layers
demlist = readRDS(paste0(path, "demlist.rds"))
chililist = readRDS(paste0(path, "chililist.rds"))
pptlist = readRDS(paste0(path, "pptlist.rds"))
tminlist = readRDS(paste0(path, "tminlist.rds"))
tmaxlist = readRDS(paste0(path, "tmaxlist.rds"))

sagelist = readRDS(paste0(path,"usgs_sagebrush/sagelist.rds"))


# --- combine environmental covariates into a summary data frame
# covariates: dem (mean,sd); chiili (mean,sd); climate (ppt, tmin, tmax) means;
dfEnv <- data.frame(dem_mean = rep(NA, nrow(fire_polygons)),
                 dem_sd = NA, chili_mean = NA, chili_sd = NA,
                 ppt = NA, tmax = NA, tmin = NA)
for(i in 1:N){
  dfEnv[i,1] = cellStats(demlist[[i]], stat = 'mean')
  dfEnv[i,2] = cellStats(demlist[[i]], stat = 'sd')
  dfEnv[i,3] = cellStats(chililist[[i]], stat = 'mean')
  dfEnv[i,4] = cellStats(chililist[[i]], stat = 'sd')
  dfEnv[i,5] = cellStats(pptlist[[i]], stat = 'mean')
  dfEnv[i,6] = cellStats(tmaxlist[[i]], stat = 'mean')
  dfEnv[i,7] = cellStats(tminlist[[i]], stat = 'mean')
}
var <- rep(NA, N)
for(i in 1:N){var[i] = cellStats(demlist[[i]], stat = "max")}
mean(var[idx])
# --- calculate pairwise mahalanobis distance
library(StatMatch)
idx <- rep(NA, N)

for(i in 1:N){idx[i] = is.finite(cellStats(sagelist[[i]][[1]], stat = 'mean'))}
n <- sum(idx)
dfEnv <- dfEnv[idx, ]

mahdist <- mahalanobis.dist(as.matrix(dfEnv))
diag(mahdist) <- NA

# --- combine spatial df with the closest matching polygons
dfspat <- fire_polygons[idx, ]
dfspat$idvar <- 1:n
rownames(mahdist) <- dfspat$idvar

dfspat$closeMatch <- NA
for(i in 1:n){dfspat$closeMatch[i] = which(mahdist[i,] == min(mahdist[i,], na.rm = T))}
dfspat$closeMatch <- apply(mahdist, 1, which.min)

#saveRDS(dfspat,paste0(path,'mtch_polygons.rds'))
# --- remove fires that are outside the 13 NLCD Landsat scenes (contain non-finite values)
#sagelist=readRDS("C:/Users/CaughlinLab/Desktop/Landsat_eros/usgs_sagebrush/sagelist.rds")

#flist=subset(sagelist,idx)
#saveRDS(flist, file=paste0(path,"usgs_sagebrush/sagelist_finite.rds"))
# --- temporarily subset sites that are small (<50)
#dfspatSubset=dfspat[dfspat$ha_clipped<50 & dfspat$Fire_Year>1990,]
#sageSubset = subset(sagebrushlist,dfspat$idvar %in% dfspatSubset$idvar)

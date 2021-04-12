library(raster)
library(sf)
library(rgdal)
library(foreach)
library(doParallel)
library(dplyr)
library(spdplyr)
################################ This script collects env data and matches sites using Mahalanobis distance metric

# --- import fire  polygons
path_fires <- "C:/Users/CaughlinLab/Desktop/Landsat_eros/nlcd_geospatial_RxFire/"

fpols <- readOGR(paste0(path_fires, list.files(path_fires, pattern = "\\_blm.shp$"))) %>%
  filter(clipAre > 10)

N <- nrow(fpols)

# --- import rasters to crop --- THIS SECTION OF CODE DOES NOT NEED TO BE REPEATED
path_rast <- "C:/Users/CaughlinLab/Desktop/Landsat_eros/nlcd_cov_raster_data/"
rlist <- list.files(path_rast, pattern = "\\.tif$")

# --- 1. Topographic data 
# Digital Elevation Model: USGS SRTM (30m): NASA/NASADEM_HGT/001
fdem <- paste0(path_rast, grep(pattern = "DEM", rlist, value = TRUE))

# CHILI: Continuous Heat-Insulation Load Index (30m): CSP/ERGo/1_0/US/CHILI
fchili <- paste0(path_rast, grep(pattern = "CHILI", rlist, value = TRUE))

# --- 2. Climate data: PRISM OREGONSTATE/PRISM/Norm81m
# Native 2.5 arcmin resolution, downsampled to 250m. 

# annual precipitation (cumulative)
fppt <- paste0(path_rast, grep(pattern = "PPT", rlist, value = TRUE))

# Minimum monthly temperature (C)
ftmin <- paste0(path_rast, grep(pattern = "TMIN", rlist, value = TRUE))

# Maximum monthly temperature (C)
ftmax <- paste0(path_rast, grep(pattern = "TMAX", rlist, value = TRUE))

# --- 3. Soil Data
# Soil Organic Carbon (250m): OpenLandMap/SOL/SOL_ORGANIC-CARBON_USDA-6A1C_M/v02
fsom <- paste0(path_rast, grep(pattern = "POLARIS", rlist, value = TRUE))

# Soil Moisture (10km): NASA_USDA/HSL/SMAP10KM_soil_moisture
# annual SD
fsmann <- paste0(path_rast, grep(pattern = "MOISTANOMSD", rlist, value = TRUE))
# March SD
fsmmar <- paste0(path_rast, grep(pattern = "MOIST03ANOMSD", rlist, value = TRUE))

# --- create a nested list for each covariate x site; first entry of a list is a name
covlist <- list()
rlist1 <- c(rlist[3], rlist[2], rlist[7])
fwgs <- spTransform(fpols, raster(paste0(path_rast, rlist1[1]))@crs)

for(i in 1:length(rlist)){
  f = raster( paste0(path_rast, rlist[i]) )
  covlist[[i]] = rlist[[i]]
  for(j in 1:N){
    covlist[[i]] = append(covlist[[i]], mask(crop(f, fwgs[j, ]), fwgs[j, ] ))
  }
}



# --- sagebrush rasters
# clip rasters to GB_LCC
# need to loop through 3 directories by years
# dcew = readOGR("C:/Users/CaughlinLab/Desktop/Landsat_eros/DryCreek/dcew_burned_area.shp")
# dcewproj <- spTransform(dcew, sagecrs)
# path_shrub <- "D:/NLCD_data_Homer_et_al_2020/sagebrush/"
# slist <- list.files(path_shrub, pattern = "\\.img$")
# 
# r <- raster(paste0(path_shrub, slist[1]))
# ts = mask(crop(r, dcewproj), dcewproj)
# ts[values(ts) > 99] = NA
# for(i in 2:length(slist)){
#   r <- raster(paste0(path_shrub, slist[i]))
#   temp = mask(crop(r, dcewproj), dcewproj)
#   temp[values(temp) > 99] = NA
#   ts <- addLayer(ts, temp)
# }
# writeRaster(ts, "C:/Users/CaughlinLab/Desktop/Landsat_eros/DryCreek/DCEW_1985_2018_sagebrush.tif")

path_sage <- "D:/NLCD_data_Homer_et_al_2020/sagebrush/"
slist <- list.files(path_sage, pattern = "\\.img$")
sagecrs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0
+ellps=WGS84 +towgs84=0,0,0,-0,-0,-0,0 +units=m +no_defs"


fproj <- spTransform(fpols, sagecrs)

# combine annual rasters
temp <- raster(paste0(path_sage, slist[1]))
for(i in 1:N) { 
  temp = mask(crop(r, fproj[i,]), fproj[i,])
  temp[values(temp) > 99] = NA
  sagelist[[i]] = temp
}
# 
for(j in 2:length(slist)) {
  r <- raster(paste0(path_sage, slist[j]))
  temp <- addLayer(temp, r)
}
sagestack <- temp

# crop and mask shagerush and shrub rasterstacks ##################################
cl = makeCluster(16)
registerDoParallel(cl)

flist = foreach(i=1:N, .packages = c('raster','rgdal', 'sf')) %dopar% {
  temp = mask(crop(sagestack, fproj[i,]), fproj[i, ])
  values(temp)[values(temp) > 99] <- NA
  temp
}
stopCluster(cl)
######################################################

# create an index for fires to be removed (no data or 0's)
ins <- rep(NA, times = N)
for(i in 1:N) {
  ins[i] = ifelse(flist[[i]][[1]]@data@max < 1 , 0, 1)
}


# ---------------------------------------------------------------- PROCEED FROM HERE
path <- "C:/Users/CaughlinLab/Desktop/Landsat_eros/"
### === Import clipped layers
# saveRDS(sagelist, paste0(path, "usgs_sagebrush/sagelist.rds"))
# sagelist = readRDS(paste0(path,"usgs_sagebrush/sagelist.rds"))


# --- combine environmental covariates into a summary data frame
dfEnv <- data.frame(dem_mean = rep(NA, N), 
                    ppt = NA, tmax = NA, tmin = NA,
                    sm03 = NA, smAnn = NA)
for(i in 1:N){
  dfEnv[i,1] = cellStats(covlist[[3]][[i+1]], stat = 'mean')
  dfEnv[i,2] = cellStats(covlist[[1]][[i+1]], stat = 'mean')
  dfEnv[i,3] = cellStats(covlist[[8]][[i+1]], stat = 'mean')
  dfEnv[i,4] = cellStats(covlist[[9]][[i+1]], stat = 'mean')
  dfEnv[i,5] = cellStats(covlist[[4]][[i+1]], stat = 'mean')
  dfEnv[i,6] = cellStats(covlist[[7]][[i+1]], stat = 'mean')
}
# site level covariate data frame

dfEnv = readRDS("dfEnv_covars.rds")

# ========================= subset fires, sage rasters, covariates based on 'outs'
n <- sum(ins)
fpols1 <- fpols[ins == T, ]
fproj1 <- fproj[ins == T, ]
dfEnv1 <- dfEnv[ins == T, ]
sagelist <- flist[ins==T]

pxldemcov <- covlist[[1]][c(0, ins) == T]
pxlchilicov <- covlist[[2]][c(0, ins) == T]
pxlsomcov <- covlist[[3]][c(0, ins) == T]

# export datasets
saveRDS(sagelist, file = paste0(path, "sagelist.rds"))
saveRDS(fproj1, file = paste0(path, "fpolygons.rds"))
saveRDS(dfEnv1, file = paste0(path, "dfEnv_covars.rds"))
saveRDS(list(dem = pxldemcov,chili =  pxlchilicov, som = pxlsomcov), file = paste0(path, "pxlcovlist.rds"))


# --- calculate pairwise Mahalanobis distance at Site Level
library(StatMatch)

mahdist <- mahalanobis.dist(as.matrix(dfEnv1))
diag(mahdist) <- NA
# Mahalanobis distance matrix
# saveRDS(mahdist, file = "mahdist_matrix.rds")


# --- combine spatial df with the closest matching polygons
# dfspat <- fire_polygons[idx, ]
# dfspat$idvar <- 1:n
# rownames(mahdist) <- dfspat$idvar
# 
# for(i in 1:n){
#   dfspat$closeMatch[i] = which(mahdist[i,] == min(mahdist[i,], na.rm = T))
# }
# dfspat$closeMatch <- apply(mahdist, 1, which.min)




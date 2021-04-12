# 
pkgs <- c("raster", "sf", "rgdal", "cluster", "spdplyr", "rgeos", "snow")
sapply(pkgs, require, character.only = T)

path <- "C:/Users/CaughlinLab/Desktop/Landsat_eros/nlcd_geospatial_RxFire/"
prjcrs <- CRS("+init=epsg:5070")

# === load base layers
gblcc <- readOGR(paste0(path, "Great Basin Landscape Conservation Cooperative Boundary/GreatBasin_LCC.shp")) %>%
  spTransform(prjcrs)

blmlcc <- readOGR(paste0(path, "BLM_National_Surface_Management_Agency/sma_wm_GreatBasin_LCC_blm/sma_wm_GreatBasin_LCC_blm.shp")) 

blmlcc@proj4string <- CRS("+init=epsg:3857")
blmlcc <- spTransform(blmlcc, prjcrs)

# === WILDFIRE DATASET
firesOne <- readOGR(paste0(path, "GB_wildfire_subset_nlcd_proj/GB_wildfire_subset/GB_wildfire_subset_1951_2018_ct_1.shp")) %>%
  filter(FireYear < 2008 & FireYear > 1986) 


firesOnebuff <- buffer(firesOne, width = -500, dissolve = FALSE) %>% spTransform(prjcrs)

f = spTransform(firesOnebuff, prjcrs)

# === LTDL DATASET

# polygons
trtmbuff <- readOGR(paste0(path, "LTDL_May_2020_Geodatabase/ltdl2020_trtm_polygons/ltdl2020_trtm_polygons_treatment_info_join_500buff.shp")) %>%
  spTransform(prjcrs)
tins <- gIntersects(trtmbuff, gblcc, byid = T) %>% as.vector()
trtmbuff <- trtmbuff[tins, ]

projbuff <- readOGR(paste0(path, "LTDL_May_2020_Geodatabase/ltdl2020_proj_polygons/ltdl2020_proj_polygons_500buff.shp")) %>%
  spTransform(prjcrs)

pins <- gIntersects(projbuff, gblcc, byid = T) %>% as.vector()
projbuff <- projbuff[pins, ]

# info
projInfo <- read.csv(paste0(path, "LTDL_May_2020_Geodatabase/ltdl_project_info.csv"))
trtmInfo <- read.csv(paste0(path, "LTDL_May_2020_Geodatabase/ltdl_trtm_info.csv"))

# points
trtmpts <- st_read(paste0(path, "LTDL_May_2020_Geodatabase/ltdl2020_trtm_pts/ltdl2020_trtm_pts.shp")) %>%
  left_join(trtmInfo, by = "Trt_ID") %>% st_transform(crs = prjcrs) %>% as_Spatial()
projpts <- st_read(paste0(path, "LTDL_May_2020_Geodatabase/ltdl2020_proj_pts/ltdl2020_proj_pts.shp")) %>%
  left_join(projInfo, by = "Prj_ID") %>% st_transform(crs = prjcrs) %>% as_Spatial()


# === Clip the fire polygons
# trtm and proj polygons

# faster in a loop with large datasets, easier to troubleshoot
fb <- firesOnebuff
for(i in 1:nrow(projbuff)) {
  print(i)
  fb = fb - projbuff[i, ]
}
fbb <- fb

for(i in 1:nrow(trtmbuff)) {
  print(i)
  fbb = fbb - trtmbuff[i, ]
}

# points
sub = fbb
sub@proj4string <- prjcrs

sub <- sub[is.na(over(sub, geometry(trtmpts))), ] # treated
sub <- sub[is.na(over(sub, geometry(projpts))), ] # projects

writeOGR(sub, paste0(path, "GBLCC_1987_2007_ct_1_Rx0.shp"), "GBLCC_1987_2007_ct_1_Rx0", driver = "ESRI Shapefile", overwrite_layer = TRUE)
# blm lands
blmlist <- list()
for(i in 1:nrow(blmlcc)) { blmlist[[i]] = buffer(blmlcc[i, ], 0) }


fires.rx0.blm <- sub
for(i in 1:length(blmlist)){
  print(i)
  fires.rx0.blm <- fires.rx0.blm - blmlist[[i]]
}

fires.rx0.blm@proj4string <- prjcrs
fires.rx0.blm$clipArea <- gArea(fires.rx0.blm, byid = TRUE) / 10000

writeOGR(fires.rx0.blm, paste0(path, "GBLCC_1987_2007_ct_1_Rx0_blm.shp"), "GBLCC_1987_2007_ct_1_Rx0_blm", driver = "ESRI Shapefile", overwrite_layer = TRUE)




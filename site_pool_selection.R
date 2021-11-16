# 
pkgs <- c("sf", "dplyr", "snow", "foreach", "doParallel")
sapply(pkgs, require, character.only = T)

path <- "D:/Landsat_eros/nlcd_geospatial_RxFire/"
prjcrs <- CRS("+init=epsg:5070")
prjepsg <- 5070

# === load base layers
gblcc <- st_read(paste0(path, "Great Basin Landscape Conservation Cooperative Boundary/GreatBasin_LCC.shp")) %>%
  st_transform(prjepsg)

blmlcc <- st_read(paste0(path, "BLM_National_Surface_Management_Agency/sma_wm_GreatBasin_LCC_blm/sma_wm_GreatBasin_LCC_blm.shp")) %>% 
  st_set_crs(3857) %>%
  st_transform(prjepsg)

# === WILDFIRE DATASET
fires <- st_read(paste0(path, "GB_wildfire_subset_nlcd_proj/GB_wildfire_subset/GB_wildfire_subset_1951_2018.shp")) %>%
  # filter(FireYear > 1986) %>% 
  filter(FireYear > 1950) %>% 
  st_transform(prjepsg) %>%
  st_buffer(0)

# https://r-spatial.github.io/sf/reference/geos_binary_ops.html 
fint <- st_intersection(fires)

mult <- fint %>% 
  filter(n.overlaps > 1) %>% 
  st_buffer(dist = 300)

multun <- st_combine(st_union(mult))


ones <- fint %>% 
  filter(n.overlaps == 1) %>%
  st_buffer(dist = 0) 

onesout <- st_difference(ones, multun)

st_write(st_collection_extract(onesout, "POLYGON"),
         paste0(path, "fire_spatial_inputs/fires_ones_1951_epsg5070.shp"), 
         delete_layer = FALSE)

# === LTDL DATASET

# polygons
trtmbuff <- st_read(paste0(path, "LTDL_May_2020_Geodatabase/ltdl2020_trtm_polygons/ltdl2020_trtm_polygons_treatment_info_join_500buff.shp")) %>%
  st_transform(prjepsg)
tins <- st_intersects(trtmbuff, gblcc)
trtmbuff <- trtmbuff[ lengths(tins) > 0, ]

projbuff <- st_read(paste0(path, "LTDL_May_2020_Geodatabase/ltdl2020_proj_polygons/ltdl2020_proj_polygons_500buff.shp")) %>%
  st_transform(prjepsg)

pins <- st_intersects(projbuff, gblcc)
projbuff <- projbuff[ lengths(pins) > 0, ]

# info
projInfo <- read.csv(paste0(path, "LTDL_May_2020_Geodatabase/ltdl_project_info.csv"))
trtmInfo <- read.csv(paste0(path, "LTDL_May_2020_Geodatabase/ltdl_trtm_info.csv"))

# points
trtmpts <- st_read(paste0(path, "LTDL_May_2020_Geodatabase/ltdl2020_trtm_pts/ltdl2020_trtm_pts.shp")) %>%
  left_join(trtmInfo, by = "Trt_ID") %>% 
  st_transform(prjepsg) 
projpts <- st_read(paste0(path, "LTDL_May_2020_Geodatabase/ltdl2020_proj_pts/ltdl2020_proj_pts.shp")) %>%
  left_join(projInfo, by = "Prj_ID") %>% 
  st_transform(prjepsg) 


# === Clip the fire polygons
# --- trtm and proj polygons
fb <- onesout

trtmbuffun <- st_combine(st_union(trtmbuff))
fbout <- st_difference(fb, trtmbuffun)

projbuffun <- st_combine(st_union(projbuff))
fbb <- st_difference(fbout, projbuffun)
# --- points
sub <- fbb

idx <- lengths(st_intersects(sub, trtmpts)) > 0 
sub <- sub[!idx, ]
idx <- lengths(st_intersects(sub, projpts)) > 0 
sub <- sub[!idx, ]

# --- write intermediate output
st_write(sub, paste0(path, "fire_spatial_inputs/fires_ones_rx0_1950_epsg5070.shp"), delete_layer = FALSE)

sub <- st_read(paste0(path, "fire_spatial_inputs/fires_ones_rx0_1950_epsg5070.shp"))

fires.rx0 <- sub

# === blm lands
blmlccun <- st_combine(st_union(blmlcc))

fires.rx0.blm <- st_intersection(fires.rx0, blmlccun)

st_write(fires.rx0.blm, paste0(path, "fire_spatial_inputs/fires_ones_rx0_blm_1950_epsg5070.shp"), 
         delete_layer = TRUE)
plot(fires.rx0.blm[,1])

findf <- fires.rx0.blm %>% 
  mutate(areaHa_clip = as.numeric(st_area(.)/1e4) ) %>% 
  filter(FireYear > 1986, FireYear < 2008) %>%
  filter(areaHa_clip > quantile(areaHa_clip, .66)) %>% 
  filter(areaHa_clip > quantile(areaHa_clip, .01), 
         areaHa_clip < quantile(areaHa_clip, .99))

# === explore in leaflet
library(leaflet)
findf %>% 
  st_transform(4326) %>% 
  leaflet() %>% addTiles() %>% 
  addPolygons()
  
# === export final product 
st_write(findf, paste0(path, "fire_spatial_inputs/GBLCC_1987_2007_ct_1_rx0_blm.shp"),
         delete_layer = TRUE)





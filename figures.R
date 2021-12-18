pckg <- c("ggplot2", "dplyr", "tidyr", "sf", "rasterVis", "ggthemes", "patchwork")
sapply(pckg, require, character.only = T)

N <- 430

# === load in data
tfires <- readRDS("data/tfires.rds")
tsage <- readRDS("data/tsage.rds")
tpxlcov <- readRDS("data/tpxlcov.rds")
tdfEnv <- readRDS("data/tdfEnv_covars.rds")
txys <- readRDS("data/txys.rds")

# === Figure 1: maps
# --- region
f <- list.files("../", pattern = "GreatBasin_LLC.shp", recursive = TRUE, full.names = TRUE)
f1 <- list.files("../", pattern = "cb_2018_us_state_5m.shp$", recursive = TRUE, full.names = TRUE)
gb <- st_read(f) %>% st_geometry()
states <- st_read(f1) %>% st_transform(4326) %>% 
  st_geometry() %>% st_crop(st_buffer(gb, 1e5))
fpts <- readRDS("data/tfires.rds") %>% 
  st_buffer(0) %>% 
  st_transform(4326) %>% 
  st_geometry() %>% 
  st_centroid() 
dat = data.frame(Area_ha = tfires$arH_clp, 
                 x = st_coordinates(fpts)[,1], 
                 y = st_coordinates(fpts)[,2])
p1 = ggplot(dat) + 
  geom_point(aes(x = x, y = y, fill = Area_ha, size = Area_ha), shape = 21, alpha = .7) +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_viridis_c(guide = "legend") +
  scale_size_continuous(range = c(0, 5)) +
  #guide_legend(colour = guide_legend(), size = guide_legend())
  geom_sf(data = gb, colour = "black", fill = NA) +
  geom_sf(data = states, colour = "grey45", fill = NA, inherit.aes = TRUE) +
  coord_sf(xlim = c(-122, -111), ylim = c(36, 45), clip = "on") +
  theme_bw() +
  theme(legend.position = "top") + #, plot.margin = margin(1,1,1,0)) + #theme(legend.key.size = unit(1, "cm")) +
  annotate("text", x = -122, y = 45, label = "(a)", size = 5)

# --- fire case #16
txys <- readRDS("data/txys.rds")
sagecrs <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0
+datum=WGS84 +units=m +no_defs"
crs_wgs84 <- st_crs(4326)

r <- tsage[[16]][,2]
xy <- txys[[16]]
rin <- rasterFromXYZ(data.frame(x = xy$x, y = xy$y, z = r), crs = sagecrs) %>% 
  projectRaster(crs = CRS(crs_wgs84$wkt)) %>% 
  as("SpatialPixelsDataFrame") %>% as.data.frame() %>% 
  rename("Cover" = "z")

ggplot() +  
  geom_tile(data=rin, aes(x=x, y=y, fill=Cover), alpha=0.8) + 
  scale_fill_viridis_c() +
  coord_equal() +
  theme_map() +
  theme(legend.position="right") +
  #theme(legend.key.width=unit(1, "cm"), legend.key.height=unit(1.5, "cm")) +
  annotate("text", x = -113.8, y = 37.975, label = "(b)", size = 5)

# --- clusters
p3 = rasterFromXYZ(data.frame(x = xy$x, y = xy$y, z = tpxlcov[[16]][,'cluster2']), crs = sagecrs) %>% 
  projectRaster(crs = CRS(crs_wgs84$wkt)) %>% 
  as("SpatialPixelsDataFrame") %>% as.data.frame() %>% 
  rename("Cluster" = "z") %>% 
  mutate(Cluster = as.factor(1)) %>% 
  ggplot() +  
  geom_tile(aes(x=x, y=y, fill=Cluster), alpha=0.8) + 
  scale_fill_viridis_d() +
  coord_equal() +
  theme_map() +
  theme(legend.position="none") +
  annotate("text", x = -113.798, y = 37.975, label = "(c)", size = 5)

p4 = rasterFromXYZ(data.frame(x = xy$x, y = xy$y, z = tpxlcov[[16]][,'cluster2']), crs = sagecrs) %>% 
  projectRaster(crs = CRS(crs_wgs84$wkt)) %>% 
  as("SpatialPixelsDataFrame") %>% as.data.frame() %>% 
  rename("Cluster" = "z") %>% 
  mutate(Cluster = as.factor(round(Cluster))) %>% 
  ggplot() +  
  geom_tile(aes(x=x, y=y, fill=Cluster), alpha=0.8) + 
  scale_fill_viridis_d() +
  coord_equal() +
  theme_map() +
  theme(legend.position="none")

p5 = rasterFromXYZ(data.frame(x = xy$x, y = xy$y, z = tpxlcov[[16]][,'cluster4']), crs = sagecrs) %>% 
  projectRaster(crs = CRS(crs_wgs84$wkt)) %>% 
  as("SpatialPixelsDataFrame") %>% as.data.frame() %>% 
  rename("Cluster" = "z") %>% 
  mutate(Cluster = as.factor(round(Cluster))) %>% 
  ggplot() +  
  geom_tile(aes(x=x, y=y, fill=Cluster), alpha=0.8) + 
  scale_fill_viridis_d() +
  coord_equal() +
  theme_map() +
  theme(legend.position="none") 
  

p6 = rasterFromXYZ(data.frame(x = xy$x, y = xy$y, z = tpxlcov[[16]][,'cluster8']), crs = sagecrs) %>% 
  projectRaster(crs = CRS(crs_wgs84$wkt)) %>% 
  as("SpatialPixelsDataFrame") %>% as.data.frame() %>% 
  rename("Cluster" = "z") %>% 
  filter(Cluster < 9) %>% 
  mutate(Cluster = as.factor(round(Cluster))) %>% 
  ggplot() +  
  geom_tile(aes(x=x, y=y, fill=Cluster), alpha=0.8) + 
  scale_fill_viridis_d() +
  coord_equal() +
  theme_map() +
  theme(legend.position="right")


p33 <- p3 | p4 | p5 | p6
  
fig1 <- (p1 | p2) / p33 + plot_annotation(theme = theme(plot.margin = margin(1,1,-1,1)))
ggsave(plot = fig1, filename = "figures/fig1.pdf", width = 240, height = 240, units = "mm")
# === Null model: requires loading outputs/null_model_lm.rds and outputs/null_model.rds files
par(mfrow = c(1, 2))
plot(density(mae), lwd = 2, xlab = "Mean absolute error, %", main = "")
plot(mae ~ jitter(tsf), pch = 19, col = rgb(0, 0, 0, .25), 
     ylab = "Mean absolute error, %", xlab = "Time since fire, years")
mtext("Null Model", side = 3, line = -2, outer = TRUE)
dev.off()
# ---
data.frame(Gompertz = mae, LinearFit = maelm) %>% 
  pivot_longer(cols = 1:2) %>% 
  rename(MAE = value, model = name) %>% 
  ggplot(aes(MAE, group = model, fill = model)) + geom_density(adjust=1.5, alpha=.5) +
  labs(y = "Density") +
  theme_bw()
# --- global models vs pre-disturbance 
plot(mae ~ kvec, pch = 19, col = rgb(.5,0,0,.75), bty = "n", xlab = "Pre-disturbance mean, %", ylab = "MAE")
points(maelm ~ kvec, pch = 19, col = rgb(0,.5,0,.75) )
abline(0, 1, lwd = 2)

# --- null MAPE and MAE 
nullout <- readRDS("outputs/null_model.rds")
tsf <- nullout$tsf
mapenull <- rep(NA, N)

for(i in 1:N) {
  y <- nullout$datnull[[i]]
  yhat <- nullout$yhatnull[[i]]
  mapenull[i] = mean( abs(y[, tsf[i]] - yhat[, tsf[i]])/y[, tsf[i]])
}
# --- null MAE
maenull <- rep(NA, N)
for(i in 1:N) {
  y <- nullout$datnull[[i]]
  yhat <- nullout$yhatnull[[i]]
  maenull[i] = mean( abs(y[, tsf[i]] - yhat[, tsf[i]]) )
}

par(mfrow = c(1, 2), mar = c(4, 4, 4, 1))
plot(density(mapenull), bty = "n", lwd = 3, xlab = "MAPE", main = "")
plot(density(maenull), bty = "n", lwd = 3, xlab = "MAE", main = "")
mtext("Null Model", side = 3, line = -2, outer = TRUE)


# === growth rate and environemtnal covariates - site level alone
# --- coefs obj from 0-cluster model
rs <- sapply(coefs, function(x) mean(x$a))
dd <- sapply(coefs, function(x) mean(x$b))

rs[rs < -5] <- NA

par(mfrow = c(2,2), mar = c(4,4,2,1))
balpha = rgb(0,0,0,.5)
plot(rs ~ tdfEnv$dem_mean, pch = 19, col = balpha, xlab = "Elevation, m", ylab = "Growth rate")
plot(rs ~ tdfEnv$tmax, pch = 19, col = balpha, xlab = "Maximum temperature, C", ylab = "Growth rate")
plot(rs ~ tdfEnv$tmin, pch = 19, col = balpha, xlab = "Minimum temperature, C", ylab = "Growth rate")
plot(rs ~ tdfEnv$ninit, pch = 19, col = balpha, xlab = "Initial population state, %", ylab = "Growth rate")




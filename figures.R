pckg <- c("brms", "raster", "ggplot2", 
          "tidyverse", "sf", "rasterVis", "ggthemes", 
          "patchwork", "grid", "ggridges", "ggspatial")
sapply(pckg, require, character.only = T)

source("helper_fns.R")
# === load in data
tfires <- readRDS("data/tfires.rds")
tsage <- readRDS("data/tsage.rds")
tpxlcov <- readRDS("data/tpxlcov.rds")
tdfEnv <- readRDS("data/tdfEnv_covars.rds")
txys <- readRDS("data/txys.rds")

N <- nrow(tfires)

# === a theme piece for reuse with other ggplots
p0 <- theme_bw() +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16))
# ==============================================

# === Figure 1: map of the study area 
# - relies on external shapefiles for a) state, b) Great Basin, c) NorthAmerica boundaries.
# - source: (a: Conservation Biology Institute: http://databasin.org/datasets/734dae972d9f4e33bb50062e0c2c5e17)
# - source: (b: www.census.gov)
# - source: (c: USGS: https://www.sciencebase.gov/catalog/item/4fb555ebe4b04cb937751db9)

# - A: region
f <- list.files("../", pattern = "GreatBasin_LLC.shp", recursive = TRUE, full.names = TRUE)
f1 <- list.files("../", pattern = "cb_2018_us_state_5m.shp$", recursive = TRUE, full.names = TRUE)
gb <- st_read(f) %>% st_transform(4326) %>% st_geometry()
states <- st_read(f1) %>% st_transform(4326) %>% 
  st_geometry() %>% 
  st_crop(st_buffer(gb, 3e4))
fpts <- readRDS("data/tfires.rds") %>% 
  st_buffer(0) %>% 
  st_transform(4326) %>% 
  st_geometry() %>% st_centroid() 

dat <- data.frame(Area_ha = tfires$arH_clp, 
                 x = st_coordinates(fpts)[,1], 
                 y = st_coordinates(fpts)[,2])
p1 <- ggplot(dat) + 
  geom_point(aes(x = x, y = y, size = Area_ha), fill = "grey45", shape = 21, alpha = .7, colour = "grey45") +
  labs(x = "Longitude", y = "Latitude") +
  scale_size_continuous("Area, ha", range = c(0, 6)) +
  geom_sf(data = gb, colour = "black", fill = NA, inherit.aes = F) +
  geom_sf(data = states, colour = "grey45", fill = NA, inherit.aes = F) +
  coord_sf(xlim = c(-122, -111), ylim = c(36, 45), clip = "on", expand = FALSE) +
  scale_x_continuous(breaks = c(-112, -114, -116, -118, -120),
                     labels = c("-112", "-114", "-116", "-118", "-120")) +
  scale_y_continuous(breaks = (c(38, 40, 42, 44)),
                     labels = c("38", "40", "42", "44")) +
  annotation_north_arrow(location = "bl") + 
  annotation_scale(location = "br", pad_x = unit(.3, "cm"), text_cex = 1.1) +
  theme_bw() +
  theme(legend.position = "right") +
  geom_rect(xmin=-113.3, xmax=-114.25, ymin=37.6, ymax=38.25, alpha = 0, colour = "black", linetype = "dotted", size = 1)

# North-America inlet
gb0 <- st_read(f) %>% st_transform(4326) %>% st_bbox() %>% st_as_sfc()
NorthAm <- st_read(list.files("../", pattern = "boundary_l_v2.shp$", recursive = TRUE, full.names = TRUE)) %>% 
  filter(COUNTRY %in% c("USA", "MEX", "CAN", "CAN USA", "MEX USA")) %>% 
  st_crop(xmin = -3000945, ymin = -3002040, xmax = 2505228, ymax = 3506052)

p1a <- ggplot() +
  geom_sf(data = NorthAm, colour = "grey65", fill = NA, inherit.aes = F) +
  geom_sf(data = gb0, alpha = 0, col = "black", size = 1.5) +
  theme_bw() + 
  theme(text = element_blank(), axis.ticks = element_blank())


fig1 <- p1 + inset_element(p1a, left = 0, right = 0.25, bottom = 0.75, top = 0.999)

ggsave(filename = "figures/fig1.pdf", fig1, width = 16, height = 16, units = "cm")
# === end Figure 1


# === Figure 2: An example of a wildfire with a clustering scheme
# fire case #16: BERYL 
txys <- readRDS("data/txys.rds")
sagecrs <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0
+datum=WGS84 +units=m +no_defs"
crs_wgs84 <- st_crs(4326)

# - A: NLCD raster
r <- tsage[[16]][,2]
xy <- txys[[16]]
rin <- rasterFromXYZ(data.frame(x = xy$x, y = xy$y, z = r), crs = sagecrs) %>% 
  projectRaster(crs = CRS(crs_wgs84$wkt)) %>% 
  as("SpatialPixelsDataFrame") %>% as.data.frame() %>% 
  rename("Cover" = "z")

p2 <- ggplot() +  
  geom_tile(data=rin, aes(x=x, y=y, fill=Cover), alpha=0.8) + 
  scale_fill_viridis_c("Cover, %", option = "D") +
  coord_equal() + 
  theme_bw() + labs(x = "Longitude", y = "Latitude") +
  theme(legend.position="right"# , panel.border = element_rect(colour = "red", size = 1.75)
        ) 

# - B-E: clusters
p3 <- rasterFromXYZ(data.frame(x = xy$x, y = xy$y, z = tpxlcov[[16]][,'cluster2']), crs = sagecrs) %>% 
  projectRaster(crs = CRS(crs_wgs84$wkt)) %>% 
  as("SpatialPixelsDataFrame") %>% 
  as.data.frame() %>% 
  rename("Cluster" = "z") %>% 
  mutate(Cluster = as.factor(1)) %>% 
  ggplot() +  
  geom_tile(aes(x=x, y=y, fill=Cluster), alpha=0.8) + 
  scale_fill_viridis_d(option = "B") +
  coord_equal() +
  theme_map() +
  theme(legend.position="none") 

p4 = rasterFromXYZ(data.frame(x = xy$x, y = xy$y, z = tpxlcov[[16]][,'cluster3']), crs = sagecrs) %>% 
  projectRaster(crs = CRS(crs_wgs84$wkt)) %>% 
  as("SpatialPixelsDataFrame") %>% as.data.frame() %>% 
  rename("Cluster" = "z") %>% 
  mutate(Cluster = as.factor(round(Cluster))) %>% 
  ggplot() +  
  geom_tile(aes(x=x, y=y, fill=Cluster), alpha=0.8) + 
  scale_fill_viridis_d(option = "B", end = 0.5) +
  coord_equal() +
  theme_map() +
  theme(legend.position="none") 

p5 = rasterFromXYZ(data.frame(x = xy$x, y = xy$y, z = tpxlcov[[16]][,'cluster5']), crs = sagecrs) %>% 
  projectRaster(crs = CRS(crs_wgs84$wkt)) %>% 
  as("SpatialPixelsDataFrame") %>% as.data.frame() %>% 
  rename("Cluster" = "z") %>% 
  mutate(Cluster = as.factor(round(Cluster))) %>% 
  ggplot() +  
  geom_tile(aes(x=x, y=y, fill=Cluster), alpha=0.8) + 
  scale_fill_viridis_d(option = "B", end = 0.75) +
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
  scale_fill_viridis_d(option = "B", end = 1) +
  coord_equal() +
  theme_map() +
  theme(legend.position="bottom", 
        # legend.text = element_text(size = 14),
        # legend.title = element_text(size = 16),
        legend.direction = "vertical")


p33 <- p3 + p4 + p5 + p6 +
  plot_layout(ncol = 2, nrow = 2) +
  plot_annotation(theme = theme(plot.margin = margin(0,0,0,0)))
  
layout <- c(
  area(t = 0, l = 1, b = 4, r = 3),
  area(t = 0, l = 4, b = 4, r = 6)
)

fig2 <- p2 + p33 + 
  plot_layout(guides = "collect") + #design = layout,
  plot_annotation(tag_levels = "A") +
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 24))
  
ggsave(filename = "figures/fig2.pdf", fig2, width = 20, height = 11, units = "cm")
# === end Figure 2


# === Figure 3:
# --- data
dat <- readRDS("outputs/null_dat.rds")
dff <- dat$dff
df <- dat$df
df0 <- dat$df0
# --- models
ms <- readRDS("outputs/null_model_brm.rds")
tempb <- ms$tempb
tempb0 <- ms$tempb0
# --- deterministic prediction
cb <- posterior_samples(tempb)
cb0 <- exp( posterior_samples(tempb0)$b_Intercept )
mb <- exp(-cb$b_Intercept / cb$b_logx)
t <- 0:99
predb <- gomp(cb$b_Intercept, cb$b_logx, cb0, t)
predb.mu <- apply(predb, 2, mean)
predb.ci <- apply(predb, 2, quantile, probs = c(0.025, 0.975))
pred <- data.frame(t = t, mu = predb.mu, mulb = predb.ci[1,], muub = predb.ci[2,])

prednoise <- sapply(predb.mu, function(x){
  rpois(nrow(cb), exp(cb$b_Intercept + cb$b_logx*log(x) + log(x)) )
})
# --- stochastic prediction
set.seed(125)
yh <- data.frame( matrix(0, nr = nrow(cb), nc = 100) )
colnames(yh) <- 0:(ncol(yh)-1)
yh[,1] <-  cb0
for(i in 2:ncol(yh)) {yh[,i] = rpois(nrow(cb), exp(cb$b_Intercept + cb$b_logx*log(yh[,i-1]) + log(yh[,i-1])) ) }
# --- end data preparation

# --- reshape and plot the results
tmp0 <- yh %>% 
  as.data.frame() %>% 
  replace(is.na(.), 0) %>% 
  mutate(id = 1:nrow(.)) %>% 
  pivot_longer(cols = 1:100) %>% 
  rename(t = name, y = value) %>% 
  mutate(t = as.numeric(t))

tmp <- tmp0 %>% 
  filter(id %in% 1:30)

yhh <- tmp0 %>% 
  dplyr::select(-id) %>% 
  group_by(t) %>% 
  summarise(across(everything(), list(sd = sd, mu = mean) )) %>% 
  ungroup() %>%
  mutate(mulb = y_mu - y_sd, muub = y_mu + y_sd) %>% 
  mutate(mulb = ifelse(mulb < 0, 0, mulb))

pnoise <- prednoise %>% 
  as.data.frame() %>% 
  replace(is.na(.), 0) %>% 
  mutate(id = 1:nrow(.)) %>% 
  pivot_longer(cols = 1:100) %>% 
  mutate(name = gsub('^.', '', name)) %>% 
  rename(t = name, y = value) %>% 
  mutate(t = as.numeric(t)) %>% dplyr::select(-id) %>% 
  group_by(t) %>% 
  summarise(across(everything(), list(sd = sd, mu = mean) )) %>% 
  ungroup() %>%
  mutate(y_sd = 1* y_sd) %>% 
  mutate(mulb = y_mu - y_sd, muub = y_mu + y_sd) %>% 
  mutate(mulb = ifelse(mulb < 0, 0, mulb))

p1 <- dff %>% 
  dplyr::select(-x, -id) %>% 
  group_by(t) %>% 
  summarise(across(everything(), list(sd = sd, mu = mean) )) %>% 
  ungroup() %>%
  mutate(y_sd = 1* y_sd) %>% 
  mutate(mulb = y_mu - y_sd, muub = y_mu + y_sd) %>% 
  mutate(mulb = ifelse(mulb < 0, 0, mulb)) %>% 
  ggplot() +
  geom_ribbon(aes(x = t, ymin = mulb, ymax = muub), alpha=0.15, fill = "maroon") +
  geom_ribbon(pnoise, mapping = aes(x = t, ymin = mulb, ymax = muub), alpha=0.5, fill = "lightgray") +
  geom_line(tmp, mapping = aes(x = t, y = y, group = id), alpha = 1, size = 0.5, colour = "skyblue") +
  geom_line(pnoise, mapping = aes(x = t, y = y_mu), size = 2, colour = "darkblue") +
  geom_line(aes(x = t, y = y_mu), size = 2, colour = "maroon") +
  labs(y = "Cover, %", x = "Time since wildfire, years") +
  xlim(0, 100) + # ylim(0, 50) +
  theme_bw() + 
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16))

ggsave(p1, filename = "figures/fig3.png", width = 180, height = 140, units = "mm")

# === end Figure 3

# === Figure 4: MAE figure
dfvis <- data.frame(insite = insample, 
                    outsite = maemat[1, ],
                    out4 = maemat[2, ], 
                    out14 = maemat[3, ]) %>% 
  mutate(null = mae.null.mu) %>% 
  pivot_longer(cols = everything()) %>% 
  rename("model" = "name", "MAE" = "value") %>% 
  mutate(model = as.factor(model)) %>% 
  mutate(model = fct_relevel(.f = model, "insite", "out14", "out4", "outsite", "null")) %>%
  mutate(model = fct_recode(model, "Null" = "null",
                            "CV: Site-level" = "outsite",
                            "In-sample: Site-level" = "insite",
                            "CV: Site-cluster 4" = "out4",
                            "CV: Site-cluster 14" = "out14"))
p1 <- dfvis %>% 
  ggplot(aes(y=model, x=MAE, height=..density..)) + 
  ggridges::geom_density_ridges(scale=1, stat="density", panel_scaling = TRUE, trim = 1) + 
  guides(fill = FALSE) + 
  theme_bw(base_size = 14) +
  labs(x = "MAE, %", y = "Model") 

p2 <- dfvis %>% 
  # mutate(model = fct_rev(model)) %>%
  group_by(model) %>% 
  summarise_all(.funs = c("mean", "sd")) %>% 
  ggplot() + 
  geom_point(aes(x = model, y = mean),colour = cols[2], size = 4, shape = 15) +
  geom_point(aes(x = model, y = sd), colour = cols[6], size = 4, shape = 17) +
  coord_flip() +
  labs(x = "Model", y = "MAE, %") +
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())

fig3 <- p1 + p2

# ggsave("figures/fig4.pdf", fig3, width = 17, height = 9, units = "cm")
# === end of Figure 3

# === Figure 4:
yerrmat <- readRDS("outputs/bias.rds")

# - calculate bias
bmat <- matrix(NA, nr = 15, nc = 29)
yerrmat <- yerrmat$yerrmat
for(i in 1:nrow(bmat)){ bmat[i,] = apply(yerrmat[i,,], 2, median, na.rm = TRUE) }

p1 <- bmat %>% 
  as.data.frame() %>% 
  rename_with(~str_sub(., start = 2), cols = everything()) %>% 
  mutate(id = 1:n(), 
         Model = fct_reorder(as.factor(c('Region', 'Site', paste0('Cluster ', 2:14))), id) ) %>% 
  pivot_longer(cols = 1:29) %>% 
  mutate(t = as.numeric(name)) %>% 
  filter(!is.na(value))  %>% 
  ggplot() +
  geom_hline(yintercept = 0, size = 1.1, colour = "black", linetype = "dashed") +
  geom_line(aes(x = t, y = value, group = id, colour = Model), size = 1.1) +
  scale_color_viridis_d(direction = -1) +
  ylab("Bias, % cover") + xlab("Time since wildfire, years") +
  theme_bw() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16)) + 
  theme(legend.position="right", legend.text = element_text(size = 14), 
      legend.title = element_text(size = 16),
      legend.direction = "vertical")
# --- end of average bias
blmat <- matrix(NA, nr = 15, nc = 29)
bumat <- matrix(NA, nr = 15, nc = 29)
for(i in 1:nrow(blmat)){ 
  blmat[i,] = apply(yerrmat[i,,], 2, quantile, prob = 0.33, na.rm = TRUE) 
  bumat[i,] = apply(yerrmat[i,,], 2, quantile, prob = 0.66, na.rm = TRUE) 
}

lumat <- blmat %>% 
  as.data.frame() %>% 
  rename_with(~str_sub(., start = 2), cols = everything()) %>% 
  mutate(id = 1:n(), 
         Model = fct_reorder(as.factor(c('Region', 'Site', paste0('Cluster ', 2:14))), id) ) %>% 
  pivot_longer(cols = 1:29) %>% 
  mutate(t = as.numeric(name)) %>% 
  filter(!is.na(value)) %>% 
  mutate(value2 = bumat %>% 
           as.data.frame() %>% 
           rename_with(~str_sub(., start = 2), cols = everything()) %>% 
           mutate(id = 1:n(), 
                  Model = fct_reorder(as.factor(c('Region', 'Site', paste0('Cluster ', 2:14))), id) ) %>% 
           pivot_longer(cols = 1:29) %>% 
           mutate(t = as.numeric(name)) %>% 
           filter(!is.na(value)) %>% pull(value))

fig4 <- p1 + 
  geom_ribbon(filter(lumat, Model == 'Cluster 14'), 
              mapping = aes(x = t, ymin = value, ymax = value2), alpha=0.15, fill = viridis::viridis(1)) + 
  geom_ribbon(filter(lumat, Model == 'Region'), 
              mapping = aes(x = t, ymin = value, ymax = value2), alpha=0.15, fill = viridis::viridis(15)[15])  

ggsave("figures/fig4.pdf", fig4, width = 19, height = 12, units = "cm")
# --- end of Figure 4

# --- Figure 5
efin <- read.csv("outputs/error_bar_N_df.csv")

pp <- efin %>% 
  dplyr::select(Null, Site, V2, V12) %>% 
  rename("Cluster 4" = V2, "Cluster 12" = V12, Region = Null) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(var = as.factor(name)) %>% 
  mutate(Model = fct_rev(fct_relevel(var, "Region", "Site", "Cluster 4", "Cluster 12"))) %>% 
  ggplot(aes(x = value, y = Model, colour = Model, fill = Model)) + 
  geom_density_ridges(scale = 3, alpha = .5) + # geom_density(alpha = .25) +
  labs(x = "Error: cover, %", y = "Density") +
  geom_vline(xintercept = 0) +
  theme_bw() +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 14), 
        strip.text = element_text(size=14))

fig4 <- p4 + pp + plot_layout(nrow = 2) + 
  plot_annotation(tag_levels = "A")

ggsave("figures/figure4.png", fig4, width = 210, height = 180, units = "mm")
# === end of Figure 4

# === Figure S4 based on Table 2: 
# summary of accuracy (out-of-sample)
# exported from prederr_sensitivity.R as .csv in /Figures/
# -- show as a figure: Figure 4
p <- read.csv("figures/Table2a.csv") %>% 
  rename(Region = Null) %>% 
  mutate(Validation = "In-sample") %>% 
  bind_rows(read.csv("figures/Table2b.csv") %>% 
              mutate(Validation = "Out-of-sample") %>% 
              rename(Region = Null)) %>% 
  pivot_longer(cols = 2:16) %>% 
  filter(name != "X2") %>% 
  mutate(name = gsub("X", "", name), var = ifelse(Metric == "MAE%" | Metric == "RMSE%", "Proportional", "Absolute")) %>%
  mutate(name = as.factor(name))  %>% 
  mutate(name = fct_relevel(name, "Region", "Site", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"))

# ===== for Figure 4
p4 <- p %>% 
  filter(var == "Absolute", Metric != "SD") %>% 
  ggplot() +
  geom_point(aes(x = name, y = value, colour = Metric, shape = Validation), size = 3, alpha = .75) +
  scale_x_discrete(labels = c("Region", "Site", paste("Cluster", c(3:14)))) +
  labs(x = "Model", y = "Error, %") +
  theme_bw() +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),) + 
  theme(legend.position="right", legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14),
        strip.text = element_text(size=14),
        legend.direction = "vertical")
# ===============

p %>% 
  filter(var == "Proportional", Metric != "SD") %>% 
  ggplot() +
  geom_point(aes(x = name, y = value, colour = Metric, shape = Validation), size = 3, alpha = .75) +
  scale_x_discrete(labels = c("Region", "Site", paste("Cluster", c(3:14)))) +
  labs(x = "Model", y = "Error, %") +
  theme_bw() +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),) + 
  theme(legend.position="right", legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14),
        strip.text = element_text(size=14),
        legend.direction = "vertical")


ggsave("figures/figureS4.pdf", plot = p, width = 210, height = 120, units = "mm")
# --- end of Figure S4

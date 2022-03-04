pckg <- c("brms", "ggplot2", "dplyr", "tidyr", "stringr", "sf", "rasterVis", "ggthemes", "patchwork", "grid")
sapply(pckg, require, character.only = T)

source("helper_fns.R")
# === load in data
tfires <- readRDS("data/tfires.rds")
tsage <- readRDS("data/tsage.rds")
tpxlcov <- readRDS("data/tpxlcov.rds")
tdfEnv <- readRDS("data/tdfEnv_covars.rds")
txys <- readRDS("data/txys.rds")

N <- nrow(tfires)

# === ggplot theme piece

p0 <- theme_bw() +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16))



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
  annotate("text", x = -121.8, y = 45, label = "(a)", size = 6)

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

p2 = ggplot() +  
  geom_tile(data=rin, aes(x=x, y=y, fill=Cover), alpha=0.8) + 
  scale_fill_viridis_c() +
  coord_equal() +
  #theme_map() +
  theme_bw() + labs(x = "Longitude", y = "Latitude") +
  theme(legend.position="right", legend.text = element_text(size = 14), 
        legend.title = element_text(size = 16)) +
  # theme(legend.key.width=unit(.5, "cm"), legend.key.height=unit(3, "cm")) +
  annotate("text", x = -113.795, y = 37.975, label = "(b)", size = 6)

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
  annotate("text", x = -113.795, y = 37.975, label = "(c)", size = 6)

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
  theme(legend.position="bottom", legend.text = element_text(size = 14), 
        legend.title = element_text(size = 16),
        legend.direction = "horizontal")


p33 <- p3 | p4 | p5 | p6 + plot_annotation(theme = theme(plot.margin = margin(0,0,0,0)))
  
fig1 <- (p1 | p2) / p33 + plot_annotation(theme = theme(plot.margin = margin(0,0,0,0))) 
  # grid.rect(width = 0.98, height = 0.98, gp = gpar(lwd = 3, col = "black", fill = NA))
ggsave(plot = fig1, filename = "figures/fig1.pdf", width = 240, height = 240, units = "mm")


# === Figure 2: null model
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
# --- end

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
  select(-id) %>% 
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
  mutate(t = as.numeric(t)) %>% select(-id) %>% 
  group_by(t) %>% 
  summarise(across(everything(), list(sd = sd, mu = mean) )) %>% 
  ungroup() %>%
  mutate(y_sd = 1* y_sd) %>% 
  mutate(mulb = y_mu - y_sd, muub = y_mu + y_sd) %>% 
  mutate(mulb = ifelse(mulb < 0, 0, mulb))

p1 <- dff %>% 
  select(-x, -id) %>% 
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

ggsave(p1, filename = "figures/fig2.pdf", width = 180, height = 140, units = "mm")

# === Figure 3: MAE figure
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

# ggsave("figures/fig3.pdf", fig3, width = 17, height = 9, units = "cm")


# === Table 2: 
# summary of accuracy (out-of-sample)
# exported from prederr_sensitivity.R as .csv in /Figures/
# -- show as a figure
p <- read.csv("figures/Table2a.csv") %>% 
  mutate(Validation = "In-sample") %>% 
  bind_rows(read.csv("figures/Table2b.csv") %>% 
              mutate(Validation = "Out-of-sample")) %>% 
  pivot_longer(cols = 2:16) %>% 
  filter(name != "X2") %>% 
  mutate(name = gsub("X", "", name), var = ifelse(Metric == "MAE%" | Metric == "RMSE%", "Proportional", "Absolute")) %>%
  mutate(name = as.factor(name))  %>% 
  mutate(name = fct_relevel(name, "Null", "Site", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13")) %>% 
  ggplot() +
  geom_point(aes(x = name, y = value, colour = Metric, shape = Validation), size = 2, alpha = .75) +
  labs(x = "Model", y = "Error") +
  facet_wrap(.~var, scales = "free", dir = "v") +
  scale_color_viridis_d(end = .8) +
  theme_bw() +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14)) + 
  theme(legend.position="right", legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14),
        strip.text = element_text(size=14),
        legend.direction = "vertical")

ggsave("figures/figure4.pdf", plot = p, width = 19, height = 12, units = "cm")
# --- end

# === Figure 3:
yerrmat <- readRDS("outputs/bias.rds")

# - calculate bias
bmat <- matrix(NA, nr = 15, nc = 29)
for(i in 1:nrow(bmat)){ bmat[i,] = apply(yerrmat[i,,], 2, median, na.rm = TRUE) }

p1 = bmat %>% 
  as.data.frame() %>% 
  rename_with(~str_sub(., start = 2), cols = everything()) %>% 
  mutate(id = 1:n(), 
         Model = fct_reorder(as.factor(c('Null', 'Site', paste0('Cluster ', 2:14))), id) ) %>% 
  pivot_longer(cols = 1:29) %>% 
  mutate(t = as.numeric(name)) %>% 
  filter(!is.na(value)) # %>% 
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
         Model = fct_reorder(as.factor(c('Null', 'Site', paste0('Cluster ', 2:14))), id) ) %>% 
  pivot_longer(cols = 1:29) %>% 
  mutate(t = as.numeric(name)) %>% 
  filter(!is.na(value)) %>% 
  mutate(value2 = bumat %>% 
           as.data.frame() %>% 
           rename_with(~str_sub(., start = 2), cols = everything()) %>% 
           mutate(id = 1:n(), 
                  Model = fct_reorder(as.factor(c('Null', 'Site', paste0('Cluster ', 2:14))), id) ) %>% 
           pivot_longer(cols = 1:29) %>% 
           mutate(t = as.numeric(name)) %>% 
           filter(!is.na(value)) %>% pull(value))

fig3 <- p1 + 
  geom_ribbon(filter(lumat, Model == 'Cluster 14'), 
              mapping = aes(x = t, ymin = value, ymax = value2), alpha=0.15, fill = viridis::viridis(1)) + 
  geom_ribbon(filter(lumat, Model == 'Null'), 
              mapping = aes(x = t, ymin = value, ymax = value2), alpha=0.15, fill = viridis::viridis(15)[15])  

ggsave("figures/fig3.pdf", fig3, width = 19, height = 12, units = "cm")
# --- end





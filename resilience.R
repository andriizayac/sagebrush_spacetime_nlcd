# ==============================================================
# This script explores the parameters in the fit models: spatial 
# and temporal variation of intrinsic growth rate and density dependence
# Note: these codes, results, and figures are not part of the forecasting manuscript
# ==============================================================


pkgs <- c("brms", "sf", "dplyr", "ggplot2", "patchwork", "lme4", "tidyr")
sapply(pkgs, require, character.only = TRUE)


# === load data
tfires <- readRDS("data/tfires.rds")
tsage <- readRDS("data/tsage.rds")
tpxlcov <- readRDS("data/tpxlcov.rds")
tdfEnv <- readRDS("data/tdfEnv_covars.rds")

# === load/transform geospatial data
omereg <- st_read("../supporting_shp_figures/us_eco_l3_state_boundaries/us_eco_l3_state_boundaries.shp") %>% 
  filter(STATE_NAME %in% c("Oregon", "Idaho",  "Wyoming",
                           "Nevada", "California", 
                           "Utah", "Arizona")) %>% 
  st_transform(4326)

fpts <- tfires %>% 
  st_buffer(0) %>% 
  st_transform(4326) %>% 
  select(which(st_agr(.) == "constant")) %>% 
  st_centroid() %>% 
  # st_geometry() %>% 
  mutate(FireNam = tfires$FireNam)

# This script shows the variation of growth rates and density-dependence across fire sites
# --- import source
modout <- readRDS("~/../../Volumes/az_drive/mae_model_outputs/modout_1.rds")
datlist <- readRDS("~/../../Volumes/az_drive/mae_model_outputs/Downloads/datlist_1.rds")

coefs <- modout$coefs
tvec <- modout$tvec

N <- length(coefs)
kvec <- sapply(tpxlcov, function(x){ mean(x[['prefire']]) })
# ---

rs <- sapply(coefs, function(x) mean(x$a) ) # growth rates
dd <- sapply(coefs, function(x) mean(x$b) ) # density-dependence
nz <- sapply(coefs, function(x) mean(x$n0) )
# rs[rs < -8] <- NA # outlier
hist(rs)
ks <- sapply(coefs, function(x) {exp(-mean(x$a)/mean(x$b))} )
ks[ks > 100] <- NA


# --- vis
p1 = data.frame(Growth.rate = rs, precover = kvec) %>% 
  ggplot() + 
  geom_point(aes(x = precover, y = Growth.rate)) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed", colour = viridis::viridis(5)[3], size = 2)+
  labs(x = "Average pre-disturbance cover, %", y = expression(paste("Growth rate [%,units ", year^{-1},"]"))) +
  theme_bw() +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16))

p2 = data.frame(Carrying.capacity = ks, precover = kvec) %>% 
  ggplot() + 
  geom_point(aes(x = precover, y = Carrying.capacity)) +
  labs(x = "Average pre-disturbance cover, %", y = "Estimated carrrying capacity, %") +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = viridis::viridis(5)[3], size = 2) + 
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16))

fig4 <- p1 + p2 
# ggsave(fig4, file = "figures/fig4.pdf", width = 25, height = 11, units = "cm")

kmae <- rep(NA, N)
for(i in 1:N) {kmae[i] = mean(abs(tpxlcov[[i]]$prefire - ks[i]))}
plot(density(kmae[!is.na(kmae)]), xlim = c(0, 20))


# === omernik stuff
# --- strange positive trend in recovery
data.frame(Growth.rate = rs, Year = tfires$FireYer) %>% 
  ggplot() +
  geom_point(aes(x = Year, y = Growth.rate)) +
  geom_abline(slope = 0, intercept = 0, size = 1.5, col = viridis::viridis(5)[3]) +
  labs(y = expression(paste("Growth rate [%,units ", year^{-1},"]"))) +
  theme_bw() +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16))

flag <- sapply(split(omereg, omereg$US_L3NAME), function(x) st_intersects(x, fpts))
flag2 <- sapply(flag, unlist)
ints <- flag2[lengths(flag2) > 0]

l3id <- names(ints)
df <- data.frame(id = 1:N, L3 = NA, gr = rs, dd = dd, ks = ks, year = tfires$FireYer)
for(i in 1:length(ints)) {
  tmp <- ints[[i]]
  df[tmp, "L3"] <- l3id[i]
}

resint <- lmer(gr ~ -1 + (1|L3), data = df)
coef(resint)
df1 <- mutate(df, year = scale(year))
restime <- lmer(gr ~ 1 + (0 + year|L3), data = df1)
coef(restime)
dfcoef <- data.frame(sl = coef(restime)$L3[,1], 
                     int = coef(restime)$L3[,2],
                     int0 = coef(resint)$L3[,1],
                     l3 = rownames(coef(restime)$L3))

osub <- omereg %>%  
  filter(US_L3NAME %in% df$L3) %>% 
  left_join(dfcoef, by = c("US_L3NAME" = "l3")) 
  
osub %>% 
  ggplot() +
  geom_sf(aes(fill = sl), colour = "black", alpha = 0.5) +
  scale_fill_viridis_c(option = "viridis", begin = 0.1) + 
  theme_bw()

# === resilience models 
# --- environment
dfres <- scale(tdfEnv) %>% 
  as.data.frame() %>% 
  mutate(tmax = scale(tmax1mu[,1])[,1]) %>% 
  bind_cols(rs = rs, dd = dd, ks = ks)


mrs <- dfres %>% 
  dplyr::select(-ks, -dd) %>% 
  brm(rs ~ ., data = ., 
      prior = prior(lasso(10), class = "b"),
      chains = 4, cores = 4, iter = 8000, warmup = 7000, seed = 123)

mdd <- dfres %>% 
  dplyr::select(-ks, -rs) %>% 
  brm(dd ~ ., data = ., 
      prior = prior(lasso(), class = "b"),
      chains = 4, cores = 4, iter = 8000, warmup = 7000)

mks <- dfres %>% 
  dplyr::select(-dd, -rs) %>% 
  brm(ks ~ ., data = ., 
      prior = prior(lasso(50), class = "b"),
      chains = 4, cores = 4, iter = 8000, warmup = 7000)



post <- posterior_samples(mrs) %>% 
  select(starts_with("b_")) %>% 
  select(-b_Intercept) %>% 
  pivot_longer(cols = everything()) %>% 
  rename(pars = name, estimate = value) %>% 
  mutate(pars = gsub("b_", "", pars)) %>% 
  mutate(pars = fct_recode(pars, "Soil Moisture Annual Anomaly: SD" = "smAnnAnomsd",
                           "Heat Load Index: mean" = "chili_mean",
                           "Elevation: SD" = "dem_sd",
                           "March Soil Moisture: mean" = "sm03mean",
                           "Heat Load Index: SD" = "chili_sd", 
                           "Annual precipitation: mean" = "ppt",
                           "March Soil Moisture: SD" = "sm03sd",
                           "Minimum temperature: mean" = "tmin", 
                           "Elevation: mean" = "dem_mean",
                           "Maximum temperature: mean" = "tmax")) %>%
  mutate(pars = fct_reorder(.f = pars, .x = estimate, .fun = mean)) %>% 
  ggplot(aes(y=pars, x=estimate, height=..density..)) + 
  ggridges::geom_density_ridges(scale=1, stat="density", panel_scaling = TRUE,trim = 1) +
  guides(fill = FALSE) + theme_bw(base_size = 14) +
  labs(x = expression(paste("Effect on population growth rate, % ", year^-1)), y = "") +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") 
  
ggsave(post, filename = "figures/fig5.pdf", width = 17, height = 9, units = "cm")

# === growth over time stan model
library(rstan)
out <- which(is.na(rs))
alist <- list(y = rs, 
              x = scale(tfires$FireYer)[,1], 
              N = length(rs),
              group = as.numeric(as.factor(df$L3)),
              k = max(as.numeric(as.factor(df$L3))) )
mvar <- stan_model("stan_code/resilience_growth_time_base.stan")
fitvar <- sampling(mvar, data = alist, chains = 4, cores = 4, iter = 8000, warmup = 7000)

summary(fitvar, pars = c("mu_a", "mu_t", "sigma_a", "sigma_t"))
plot(fitvar, pars = c("mu_a", "mu_t", "sigma_a", "sigma_t"))

postvar <- sapply(extract(fitvar, pars = c("mu_a", "mu_t", "sigma_a", "sigma_t")), rbind)

fig6c <- postvar %>% 
  as.data.frame() %>% 
  mutate(sigma_a = exp(sigma_a), sigma_t = exp(sigma_a + sigma_t)) %>% 
  pivot_longer(cols = everything()) %>% 
  rename(pars = name, estimate = value) %>%
  mutate(pars = as.factor(pars)) %>% 
  mutate(pars = fct_recode(pars, "Average Growth Rate" = "mu_a",
                           "Year effect on mean" = "mu_t",
                           "Average variance" = "sigma_a",
                           "Year effect on variance" = "sigma_t")) %>%
  ggplot(aes(y=pars, x=estimate, height=..density..)) + 
  ggridges::geom_density_ridges(scale=1, stat="density", panel_scaling = TRUE,trim = 1) +
  guides(fill = FALSE) + theme_bw(base_size = 14) +
  labs(x = "Effect size", y = "") +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  annotate("text", x = .3, y = 4.3, label = "(c)", size = 6)


fig6a <- df %>% 
  dplyr::filter(gr > -5) %>% 
  ggplot(aes(x = year, y = gr)) + geom_point(alpha = .75) +
  guides(fill = FALSE) + theme_bw(base_size = 14) +
  labs(x = "Year", y = expression(paste("Growth rate, % ", year^-1)) ) +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  annotate("text", x = 1987, y = 2, label = "(a)", size = 6)


fig6b <- bind_cols(st_coordinates(fpts), gr = df$gr) %>% 
  filter(gr > -5) %>% 
  ggplot() +
  geom_sf(data = osub, aes(fill = int0), colour = "black", alpha = 0.9) +
  scale_fill_viridis_c(option = "viridis", begin = 0.1) +
  #geom_point(aes(x = X, y = Y, size = gr), alpha = .6) +
  #scale_size_continuous(range = c(0, 2))
  coord_sf(xlim = c(-123.2, -109), ylim = c(36, 45), clip = "on") +
  theme_bw() +
  theme(legend.position = "right", axis.text = element_text(size = 7), axis.title = element_text(size =12), legend.title = element_blank()) + 
  annotate("text", x = -122.8, y = 44.65, label = "(b)", size = 6)

  

fig6a / (fig6b + fig6c)
ggsave("figures/fig6.pdf", width = 24, height = 14, units = "cm")


# === idaho level
# idahoco <- st_read("../supporting_shp_figures/tl_2016_16_cousub/tl_2016_16_cousub.shp") %>% 
#   st_transform(4326)
# 
# flagi <- sapply(split(idahoco, idahoco$NAME), function(x) st_intersects(x, fpts))
# flagi2 <- sapply(flagi, unlist)
# ints <- flagi2[lengths(flagi2) > 0]
# 
# l3id <- names(ints)
# df <- data.frame(id = 1:N, L3 = NA, gr = rs, dd = dd, ks = ks, year = tfires$FireYer)
# for(i in 1:length(ints)) {
#   tmp <- ints[[i]]
#   df[tmp, "L3"] <- l3id[i]
# }
# 
# 
# resint <- lmer(gr ~ -1 + (1|L3), data = df)
# coef(resint)
# df1 <- mutate(df, year = scale(year))
# restime <- lmer(gr ~ 1 + (0 + year|L3), data = df1)
# coef(restime)
# dfcoef <- data.frame(sl = coef(restime)$L3[,1], 
#                      int = coef(restime)$L3[,2],
#                      int0 = coef(resint)$L3[,1],
#                      l3 = rownames(coef(restime)$L3))
# 
# osub <- idahoco %>%  
#   filter(NAME %in% df$L3) %>% 
#   left_join(dfcoef, by = c("NAME" = "l3")) 
# 
# 
# ggplot() +
#   geom_sf(data = osub, aes(fill = int0), colour = "black") +
#   geom_sf(data = fpts) +
#   coord_sf(xlim = c(-117, -111.5), ylim = c(42, 45), clip = "on") +
#   scale_fill_viridis_c(option = "D",begin = 0.1) +
#   theme_bw()
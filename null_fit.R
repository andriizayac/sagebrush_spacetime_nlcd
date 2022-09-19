pkgs <- c("dplyr", "tidyr", "lme4", "ggplot2", "brms")
sapply(pkgs, require, character.only = T)

years <- c(1985:2018)

source("helper_fns.R")

# ====================================
# inputs from pxlmatching.R
tfires <- readRDS("data/tfires.rds")  # fire polygons and metadata
tsage <- readRDS("data/tsage.rds")  # [pixel,year] NLCD data
tpxlcov <- readRDS("data/tpxlcov.rds") # pixel-level covariates
tdfEnv <- readRDS("data/tdfEnv_covars.rds") # site-level covariates 
txys <- readRDS("data/txys.rds") # pixel coordinates 

N <- nrow(tfires)

# === create a directory to store model outputs
dir.create("outputs/", showWarnings = FALSE)

# === create data inputs
# - select 5% of the observations from each sample 
ssize <- sapply(tsage, function(x) floor(nrow(x)*.05) )
set.seed(125)
nullsage <- lapply(tsage, function(x) {
  ind = sample(1:nrow(x), floor(nrow(x)*.1), replace = F)
  x[ind,]
})

# - initiate a combine data sets from the list
dff <- dat.gen.null.full(nullsage, tfires, 1)
df <- dat.gen.null(nullsage, tfires, 1)
df0 <- dat.gen.null.init(nullsage, tfires, 1)
for(i in 2:N) {
  dff <- rbind(dff, dat.gen.null.full(nullsage, tfires, i))
  df <- rbind(df, dat.gen.null(nullsage, tfires, i))
  df0 <- rbind(df0, dat.gen.null.init(nullsage, tfires, i))
}
# saveRDS(list(df = df, df0 = df0, dff = dff), file = "outputs/null_dat.rds")

# === brms
######################################################
nulldat <- readRDS("outputs/null_dat.rds")
ymus <- list()
for(i in 1:N) {
  ymus[[i]] <- tsage[[i]][, c(tfires$FireYer[i]-1984):31] %>% 
    summarize_all(mean) %>% 
    as.matrix()
}
# === using Bayesian estimation to get pars uncertainty
library(brms)
input <- readRDS("outputs/null_dat.rds")
bdff <- input$dff
bdf <- input$df
bdf0 <- input$df0
# --- fit brms models
# tempb <- brm(y ~ 1 + log(x) + offset(log(x)),
#            data = bdf,
#            family = poisson,
#            algorithm = "sampling", chains = 1, cores = 1, iter = 100)
# 
# tempb0 <- brm(cover ~ 1 + t,
#              data = bdat0,
#              family = poisson,
#              algorithm = "sampling", chains = 1, cores = 1)
# ---
ms <- readRDS("outputs/null_model_brm.rds")
tempb <- ms$tempb
tempb0 <- ms$tempb0
# ---
cb <- posterior_samples(tempb)
cb0 <- exp( posterior_samples(tempb0)$b_Intercept )
mb <- exp(-cb$b_Intercept / cb$b_logx)
t <- 0:99

predb <- gomp(cb$b_Intercept, cb$b_logx, cb0, t)
predb.mu <- apply(predb, 2, mean)
predb.ci <- apply(predb, 2, quantile, probs = c(0.025, 0.975))
matplot(t(predb), type = "l")
pred <- data.frame(yhat = predb.mu, t = t)
# === end of brms fit

# stochastic Poisson predictions
set.seed(125)
yh <- data.frame( matrix(0, nr = nrow(cb), nc = 100) )
yh[,1] <- rpois(nrow(cb), cb0)
for(i in 2:ncol(yh)) {yh[,i] = rpois(nrow(cb), exp(cb$b_Intercept + cb$b_logx*log(yh[,i-1]) + log(yh[,i-1])) ) }

matplot(t(yh), type = "l")

# === leave-one-out lme4 model predictions
# for(i in 1:N){
#   print(i)
#   dfloo <- filter(df, id != i)
#   temp <- glm(y ~ 1 + log(x), offset = log(x), family = poisson(link = "log"), data = dfloo) #
#   # --- estimate initial population size based the first 5 years of data post-fire
#   df0loo <- filter(df0, id != i)
#   temp_n0 <- glm(cover ~ 1 + t, family = poisson(link = "log"), data = df0loo) # 
# 
#   coefsloo[[i]] = list(coef(temp), coef(temp_n0))
# }

coefsloo <- readRDS("outputs/null_loo_coefs.rds")

rloo <- sapply(coefsloo, function(x) {x[[1]][[1]] } )
bloo <- sapply(coefsloo, function(x) {x[[1]][[2]] } )
n0loo <- sapply(coefsloo, function(x) {x[[1]][[1]] } )

# dfyhat = data.frame(yhat = gomp(rloo[1], bloo[1], exp(n0loo[1]), 0:max(bdff$t[bdff$id == 1]))[1,], id = 1, t = 0:max(bdff$t[bdff$id == 1]))
# dfyhat = data.frame(yhat = gomp(mean(cb[,1]), mean(cb[,2]), mean(exp(cb0)), 0:max(bdff$t[bdff$id == 1]))[1,], id = 1, t = 0:max(bdff$t[bdff$id == 1]))
dfyhat = data.frame(yhat = gomp(rloo[1], bloo[1], exp(n0loo[1]), 0:28)[1,], id = 1, t = 0:28)

for(i in 2:N){
   # tmp = data.frame(yhat = gomp(rloo[i], bloo[i], exp(n0loo[i]), 0:max(bdff$t[bdff$id == i]))[1,], id = i, t = 0:max(bdff$t[bdff$id == i]))
   # tmp = data.frame(yhat = gomp(mean(cb[,1]), mean(cb[,2]), mean((cb0)), 0:max(bdff$t[bdff$id == i]))[1,], id = i, t = 0:max(bdff$t[bdff$id == i]))
   tmp = data.frame(yhat = gomp(rloo[i], bloo[i], exp(n0loo[i]), 0:28)[1,], id = i, t = 0:28)
   dfyhat = bind_rows(dfyhat, tmp)
}

dfyhat %>% 
  group_by(id) %>% 
  mutate(myr = max(t)) %>% 
  ungroup() %>% 
  ggplot() + geom_line(aes(x = t, y = yhat, group = id, colour = myr), alpha = 0.5, size = 1.5) + 
  scale_color_viridis_c() +
  facet_wrap(.~myr)+
  theme_bw() + theme(legend.position = "none")

write.csv(dfyhat, file = "outputs/null_dat_pred_in.csv", row.names = FALSE)

# === calculate site-level average and match predictions for year 10
dff.fin.t9 <- bdff %>% 
  mutate(t = t + 1) %>% 
  select(-y) %>% 
  group_by(id, t) %>% 
  summarize_all(mean) %>% 
  ungroup() %>% 
  rename(y = x) %>% 
  filter(t == 8) %>%
  mutate(yhat = dfyhat %>% 
           filter(t == 9) %>% 
           pull(yhat))
write.csv(dff.fin.t9, file = "outputs/null_dat_pred_T9_in.csv")
# --- end



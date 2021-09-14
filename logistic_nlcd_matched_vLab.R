# === libraries + paths
pkgs <- c("brms", "matrixcalc", "tidyverse", "lme4", "raster", "rgdal", "spdplyr", "doParallel", "parallel", "snow")
sapply(pkgs, require, character.only = T)

years <- c(1985:2018)

# ====================================
# see pxlmatching.R 
tfires <- readRDS("data/tfires.rds")
tsage <- readRDS("data/tsage.rds")
tpxlcov <- readRDS("data/tpxlcov.rds")
tdfEnv <- readRDS("data/tdfEnv_covars.rds")

kvec <- sapply(tpxlcov, function(x){ mean(x[['prefire']]) })

# f <- function(x) {
#   matplot(t(x), type = "l", col = rgb(.5,0,.5,.25), main = i);
#   i <<- i+1
# }
# f(tsage[[i]])
# --- parametrize k prior
# remotes::install_github("aursiber/fitdistrplus")
library(fitdistrplus)
klist <- list()
for(i in 1:length(tsage)) { klist[[i]] = tpxlcov[[i]]$prefire }
kprior <- unlist(klist)
fgamma <- fitdist( kprior[kprior > 0], "gamma" )
# plot(density(kprior), lwd = 2)
# curve(dgamma(x, shape = fgamma$estimate[1], rate = fgamma$estimate[2]), add = T, col ="brown", lwd = 3)

# --- Poisson GLM on the log-log scale
# generate sample model for compilation
glm.dat <- function(tsage, tfires, tpxlcov, i=NULL, clN = NULL) {
  # this function generates a data input for log-log gompertz Poisson model
  # output: response - sage % for t = 2,3,..T
  #   predictor - negative sage %, encoding 0 -> .25, for t = 1,2,...T-1
  #   group - cl - pixel-level classes for random effect
  mat <- data.matrix(tsage[[i]][, c(tfires$FireYer[i]-1984):33])
  xpr <- vec(mat[,-ncol(mat)])
  dat <- data.frame(y = vec(mat[,-1]), 
                  x = -log(ifelse(xpr == 0, .25, xpr)),
                  cl = as.factor(rep(tpxlcov[[i]][,paste0("cluster", clN)], times = ncol(mat)-1)),
                  t = sort( rep(2:ncol(mat), nrow(mat)) ) 
                  )
}

dat <- glm.dat(tsage, tfires, tpxlcov, 4)
bprior <- c(prior(normal(0, .75), class = b, coef = Intercept),
            set_prior("lognormal(-2, 1.25)", lb = 0, class = "b"))

temp <- brm(y ~ 0 + Intercept + x + offset(-x), family = "poisson", prior = bprior, #sample_prior = "only",
            data=dat, chains = 4, cores = 4, iter = 100)
stanfit <- temp$fit

plot( conditional_effects(temp, method = "predict"), points = F)
pred <- predict(temp)
plot(dat$y ~ (-1*dat$x))
points(pred[,1] ~ -dat$x, pch = 19, col = rgb(0,0,0,.1), cex = .2)

# === glmer predictions
cl <- makeCluster(8)
registerDoParallel(cl)
foreach(i = c(46, 98), #100:length(tsage), 
                     .packages = c('matrixcalc', 'brms'),
                     .export = c("tsage", "tfires", "tpxlcov", "stanfit", "path", "glm.dat")) %dopar% {
  # extract data
  dat <- glm.dat(tsage, tfires, tpxlcov, i)
  # prior
  bprior <- c(prior(normal(0, .85), class = b, coef = Intercept),
              set_prior("lognormal(-2, 1.25)", lb = 0, class = "b"))
  # fit brms model
  temp <- brm(y ~ 0 + Intercept + x + offset(-x), family = "poisson", prior = bprior, data = dat,
              chains = 4, cores = 4,  iter = 2000, warmup = 1800, control = list(adapt_delta = 0.95))
  saveRDS(temp, file = paste0(path, "models_stan_eros/gomp_glm_fe_", i, ".rds"))
  rm(temp)
}
stopCluster(cl)
# in-sample predictions
modlist <- list.files(paste0(path, "/models_stan_eros/"), pattern = "\\.rds")
modlist <- modlist[order(nchar(modlist), modlist)]
cl <- makeCluster(23)
registerDoParallel(cl)
predlist <- foreach(i = 1:length(modlist), .packages = c('brms'), .export = c("modlist")) %dopar% {
          m <- readRDS(paste0(path, "/models_stan_eros/", modlist[i]))
          return(predict(m))
        }
stopCluster(cl)

# === lmer models and predictions ####
# === glmer 
mlist.glm <- list()
mlist0.glm <- list()
coefs <- list()
datlist.glm <- list()
tvec <- rep(NA, length(kvec))

for(i in 1:3){ #length(tsage)
  print(i)
  
  # estimate growth and K parameters
  dat <- glm.dat(tsage, tfires, tpxlcov, i, 3)
  temp <- lmer(y ~ (1|cl) + (0+x|cl), REML = TRUE, data = dat) #
  
  # estimate initial population size based the first 5 years of data post-fire
  dat.init <- glm.dat.init(tsage, tfires, tpxlcov, i, 3)
  temp_n0 <- glmer(cover ~ (1|cl) + (0+t|cl), family = Gamma(link = "log"), data = dat.init) # 
  
  # store models
  mlist.glm[[i]] = temp
  mlist0.glm[[i]] = temp_n0
  
  datlist.glm[[i]] <- data.matrix(tsage[[i]][,c(tfires$FireYer[i]-1984):33])
  tvec[i] <- ncol( datlist.glm[[i]])
  coefs[[i]] = data.frame(a = coef(temp)$cl[,2], b = coef(temp)$cl[,1], n0 = exp(coef(temp_n0)$cl[,2]))
}

sapply(coefs, function(x) { sum(x$n0 < 0) } ) # check for negative N0

#  === Gomperz equation
# K * exp(C * exp(-alpha*t)) 
gomp <- function(a, b, n0, t) {
  # vector-valued Gomperz model
  # in: coefs from log-lin and N0-Gamma models and a vector of time steps
  # out: matrix of values nrow = length(a), and ncol = length(t)
  K <- exp(-a/b)
  C <- log(n0/K)
  u <- K * exp(C * exp( outer(b, t) ))
  return(u)
}

# === In-Sample predictions 
##################################
errvec.glm <- rep(NA, times = length(tsage))

for(j in 1:length(tsage)) {
  a <- coefs[[j]]$a
  b <- coefs[[j]]$b
  n0 <- coefs[[j]]$n0
  t <- 1:tvec[j]-1
  matpred <- gomp(a, b, n0, t) %>% 
    as.data.frame() %>% 
    mutate(cl = 1:n()) 
  
  a <- as.data.frame(datlist.glm[[j]]) %>%
    mutate(cl = tpxlcov[[j]]$cluster) %>% 
    left_join(matpred, by = c("cl" = "cl"))
  
  obs <- as.matrix( a[, 1:length(t)] )
  pred <- as.matrix( a[, c(length(t)+2):ncol(a)] )
  errvec.glm[j] <- mean( abs(obs - pred))
}

errvec.glm[is.infinite(errvec.glm) | errvec.glm > 100] <- NA
plot(errvec.glm~kvec, xlab = "Pre-fire cover", ylab = "MAE", pch = 19, main = "In-sample error")
abline(0, 1, lty = "dashed", lwd = 2, col = "gray")
##################################

# ===== Forecasting accuracy

# === Matching at site and pixel levels (pixel-to-cluster)
##################################
mah <- StatMatch::mahalanobis.dist(tdfEnv[,])
diag(mah) <- NA
match <- apply(mah, 1, which.min)
mahdist <- apply(mah, 1, min, na.rm = T)

pxl.dist <- function(dat, dat.m) {
  # finds a match for each focal pixel to a cluster in the reference site
  # - replaces cluster value for matched cluster
  # dat: focal pixel-level covariates
  # dat.m: reference pixel-level covariates
  df <- dat.m %>%
    group_by(cluster) %>%
    summarize_all(mean) %>% dplyr::select(-cluster)
  mah.p <- StatMatch::mahalanobis.dist(dat[, 1:5], df)
  
  dat$cluster <- apply(mah.p, 1, which.min)
  return(dat)
}

errvec.glm.out <- rep(NA, length(tsage))
errvec.glm.out.rand <- rep(NA, length(tsage))

err.yr.out <- matrix(NA, nr = length(tsage), nc = 31)

for(i in 1:length(tsage)){
  # glmer out-of-sample Mahalanobis nearest match
  j <- match[i]
  a <- coefs[[j]]$a
  b <- coefs[[j]]$b
  n0 <- coefs[[j]]$n0
  t <- 1:tvec[i]-1
  
  matpred <- as.data.frame(gomp(a, b, n0, t)) %>% 
    mutate(cl = row_number())
  pred.m <- pxl.dist(tpxlcov[[i]], tpxlcov[[j]]) %>% 
    dplyr::select(cluster) %>% 
    left_join(matpred, by = c("cluster" = "cl")) %>% 
    dplyr::select(-cluster)
  
  # glmer out-of-sample, random match
  l <- sample(1:length(tsage), 1, replace = T)
  a <- coefs[[l]]$a
  b <- coefs[[l]]$b
  n0 <- coefs[[l]]$n0
  t <- 1:tvec[i]-1
  
  matpred.r <- as.data.frame(gomp(a, b, n0, t)) %>% 
    slice(sample(1:n())) %>%
    mutate(cl = row_number())
  pred.r <- tpxlcov[[i]] %>%
    dplyr::select(cluster) %>% 
    left_join(matpred.r, by = c("cluster" = "cl")) %>%
    dplyr::select(-cluster)
  
  obs <- datlist.glm[[i]] 
  
  errvec.glm.out[i] <- mean( as.matrix( abs(pred.m - obs)) )
  errvec.glm.out.rand[i] <- mean( as.matrix( abs(pred.r - obs)) )
  for(k in 1:length(t)) {
    err.yr.out[i, k] = mean( (obs[, k] - pred.m[, k])/pred.m[, k], na.rm = T)
  }
}
errvec.glm.out[!is.finite(errvec.glm.out) | errvec.glm.out > 100] <- NA
errvec.glm.out.rand[!is.finite(errvec.glm.out.rand) | errvec.glm.out.rand > 100] <- NA
err.yr.out[!is.finite(err.yr.out) | err.yr.out > 100 | err.yr.out < 100] <- NA
boxplot(err.yr.out)
plot(apply(err.yr.out, 2, mean, na.rm = T))

plot(errvec.glm.out ~ mahdist, col = "purple", pch = 19, ylim = c(0, 17), 
     xlab = "Mahalanobis distance", ylab = "Mean Absolute Error, %", 
     main = "Out-of-sample accuracy")
points(errvec.glm.out.rand ~ mahdist, col = "gold", pch = 19, ylim = c(0, 12))


plot(density(errvec.glm.out, na.rm = T), lwd = 3, col = "purple", xlab = "Mean absolute error, %")
lines(density(errvec.glm.out.rand, na.rm = T), lwd = 3, col = "gold")

cat("In-sample MAE (%): ", mean(errvec.glm, na.rm = TRUE))
cat("Out-of-sample MAE (%): ", mean(errvec.glm.out, na.rm = TRUE))
cat("NULL out-of-sample MAE (%): ", mean(errvec.glm.out.rand, na.rm = TRUE))

data.frame(InSample = errvec.glm, 
           OutSample = errvec.glm.out, 
           OutSample_Null = errvec.glm.out.rand) %>% 
  pivot_longer(cols = 1:3) %>% 
  filter(value < 100) %>% 
  ggplot(aes(x = name, y = value)) + geom_boxplot() +
  labs(x = "Model matching", y = "Mean absolute error, %") + theme_classic(base_size = 14)

data.frame(kvec = kvec, mae = errvec.glm.out) %>%
  ggplot(aes(x = kvec, y = mae)) + geom_point() +
  geom_abline(slope = 1, intercept = 0, colour = "gray", size = 1.5, linetype = "dashed",  show.legend = F) +
  labs(x = "Pre-disturbance average cover, %", y = "Mean absolute error, %") +
  theme_classic(base_size = 14)


e <- as.data.frame(err.yr.out)
names(e) <- 1:ncol(e)
e %>%
  pivot_longer(1:ncol(e)) %>%
  filter(is.finite(value), !is.na(value)) %>%
  group_by(name) %>%
  mutate(mean = mean(value), upper = mean(value) + sd(value), lower = mean-sd(value)) %>%
  mutate(name = as.numeric(name)) %>%
  ungroup() %>%
  distinct(mean, .keep_all = TRUE) %>% 
  ggplot(aes(x=name, y = mean)) + geom_point(size = 3) + 
  geom_errorbar(aes(ymax = upper, ymin = lower), width = .05, size = 1) +
  geom_hline(yintercept = 0, colour = "purple", size = 1.5, linetype = "dashed",  show.legend = F) +
  labs(x = "Time since fire, years", y = expression(paste("Error (Y - ", hat(y),"), %"))) +
  theme_classic(base_size = 14)


# bind_cols(tdfEnv, MAE = errvec.glm.out) %>% GGally::ggpairs()

ggsave("FigureC.pdf", width = 150, height = 120, units = "mm", path = paste0(path, "figs/"))
##################################

# === matching at the site level only
##################################
mah <- StatMatch::mahalanobis.dist(tdfEnv[,])
diag(mah) <- NA
match <- apply(mah, 1, which.min)
mahdist <- apply(mah, 1, min, na.rm = T)


errvec.glm.out <- rep(NA, length(tsage))
errvec.glm.out.rand <- rep(NA, length(tsage))

for(i in 1:length(tsage)){
  obs <- datlist.glm[[i]]
  # glmer out-of-sample Mahalanobis nearest match
  j <- match[i]
  a <- mean(coefs[[j]]$a)
  b <- mean(coefs[[j]]$b)
  n0 <- mean(coefs[[j]]$n0)
  t <- 1:tvec[i]-1
  
  matpred.m <- matrix(gomp(a, b, n0, t), nr = nrow(obs), nc = length(t), byrow = T)
  
  # glmer out-of-sample, random match
  l <- sample(1:length(tsage), 1, replace = T)
  a <- mean(coefs[[l]]$a)
  b <- mean(coefs[[l]]$b)
  n0 <- mean(coefs[[l]]$n0)
  t <- 1:tvec[i]-1
  
  matpred.r <- matrix(gomp(a, b, n0, t), nr = nrow(obs), nc = length(t), byrow = T)
  
  errvec.glm.out[i] <- mean( abs(matpred.m - obs)) 
  errvec.glm.out.rand[i] <- mean( abs(matpred.r - obs)) 
}
errvec.glm.out[!is.finite(errvec.glm.out) | errvec.glm.out > 100] <- NA
errvec.glm.out.rand[!is.finite(errvec.glm.out.rand) | errvec.glm.out.rand > 100] <- NA

plot(errvec.glm.out ~ mahdist, col = "purple", pch = 19, ylim = c(0, 17), 
     xlab = "Mahalanobis distance", ylab = "Mean Absolute Error, %", 
     main = "Out-of-sample accuracy")
points(errvec.glm.out.rand ~ mahdist, col = "gold", pch = 19, ylim = c(0, 12))


plot(density(errvec.glm.out, na.rm = T), lwd = 3, col = "purple", xlab = "Mean absolute error, %")
lines(density(errvec.glm.out.rand, na.rm = T), lwd = 3, col = "gold")

cat("In-sample MAE (%): ", mean(errvec.glm, na.rm = TRUE))
cat("Out-of-sample MAE (%): ", mean(errvec.glm.out, na.rm = TRUE))
cat("NULL out-of-sample MAE (%): ", mean(errvec.glm.out.rand, na.rm = TRUE))

data.frame(InSample = errvec.glm, 
           OutSample = errvec.glm.out, 
           OutSample_Null = errvec.glm.out.rand) %>% 
  pivot_longer(cols = 1:3) %>% 
  filter(value < 100) %>% 
  ggplot(aes(x = name, y = value)) + geom_boxplot() +
  labs(x = "Model matching", y = "Mean absolute error, %") + theme_classic(base_size = 14)

##################################

# === multiple matches averaged
##################################

mah <- StatMatch::mahalanobis.dist(tdfEnv[,])
diag(mah) <- NA
# match is a set of models falling within 1 Std Dev from the min (i.e., nearest match)
match <- apply(mah, 1, function(x) { which(x < min(x, na.rm = T) + sd(x, na.rm = T)) } )

errvec.glm.out <- rep(NA, length(tsage))
errvec.glm.out.rand <- rep(NA, length(tsage))

for(i in 1:length(tsage)){
  a <- datlist.glm[[i]] # reference data
  
  # glmer out-of-sample Mahalanobis caliper interval mean
  j <- match[[i]]
  amean <- mean( sapply(coefs[j], function(x) { x$a }) )
  bmean <- mean( sapply(coefs[j], function(x) { x$b }) )
  n0mean <- mean( sapply(coefs[j], function(x) { x$n0 }) )
  t <- 1:tvec[i]-1
  matpred <- matrix(gomp(amean, bmean, n0mean, t)[1,], nr = nrow(a), nc = length(t), byrow = T)
  
  # glmer out-of-sample, random match
  l <- sample(1:length(tsage), length(j), replace = T)
  amean <- mean( sapply(coefs[l], function(x) { x$a }) )
  bmean <- mean( sapply(coefs[l], function(x) { x$b }) )
  n0mean <- mean( sapply(coefs[l], function(x) { x$n0 }) )

  matpred.rand <- matrix(gomp(amean, bmean, n0mean, t)[1,], nr = nrow(a), nc = length(t), byrow = T)
  
  errvec.glm.out[i] <- mean(abs(matpred - a))
  errvec.glm.out.rand[i] <- mean(abs(matpred.rand - a))
}
errvec.glm.out[!is.finite(errvec.glm.out) | errvec.glm.out > 100 ] <- NA
errvec.glm.out.rand[!is.finite(errvec.glm.out.rand) | errvec.glm.out > 100] <- NA

plot(errvec.glm.out ~ mahdist, col = "purple", pch = 19, ylim = c(0, 12))
points(errvec.glm.out.rand ~ mahdist, col = "gold", pch = 19, ylim = c(0, 12))


plot(density(errvec.glm.out, na.rm = T), col = rgb(1,.5,1,1), xlab = "MAE")
lines(density(errvec.glm.out.rand, na.rm = T),col = rgb(0,.5,.1,1))

cat("In-sample MAE (%): ", mean(errvec.glm, na.rm = T))
cat("Out-of-sample MAE (%): ", mean(errvec.glm.out, na.rm = T))
cat("NULL out-of-sample MAE (%): ", mean(errvec.glm.out.rand, na.rm = T))
##################################

######################
# === draft code
{
# --- in-sample error estimates
err.yr.glm <- matrix(NA, length(modlist), 33)
resid.yr.glm <- matrix(NA, length(modlist), 33)
err.glm <- rep(NA, length(modlist))
datlist <- list()
# calculate accuracy
for(i in 1:length(modlist)) {
  # k <- as.numeric( regmatches(modlist[i], regexec('gomp_glm_fe_(.*?)\\.', modlist[i]))[[1]][2] )
  dat <- glm.dat(tsage, tfires, tpxlcov, i)
  err.glm[i] = mean(abs(predlist[[i]][,1] - dat$y))
  for(t in 2:33) {
    # mean absolute error
    err.yr.glm[i,t] = mean(abs(predlist[[i]][dat$t == t,1] - dat$y[dat$t == t]))
    # average residuals by year
    resid.yr.glm[i,t] <- mean(predlist[[i]][dat$t == t,1] - dat$y[dat$t == t]) 
  }
  datlist[[i]] <- dat
}
plot(err.glm~kvec)
par(mfrow = c(1, 2))
boxplot(err.yr.glm[ ,1:30], main = "GLM", ylab = "MAE (% cover)", xlab = "Time since fire (years)")

boxplot(resid.yr.glm[ ,1:30], main = "GLM", ylab = "Residuals [ y_hat - y ]", xlab = "Time since fire (years)")
abline(h = .5, col = "red", lty = "dashed")

# --- out-of-sample errors; match sites at the fire level
tdfEnv$u0 <- sapply(datlist, function(x) { mean(x$y[x$t == 2]) } )

mah <- StatMatch::mahalanobis.dist(tdfEnv[,-5])
diag(mah) <- NA
match <- apply(mah, 1, which.min)
mahdist <- apply(mah, 1, min, na.rm = T)


cl <- makeCluster(20)
registerDoParallel(cl)
lout <- foreach(i = 1:length(modlist), .packages = c('brms')) %dopar% {
  # glmer out-of-sample predictions
  m.glm <- readRDS(paste0(path, "models_stan_eros/", modlist[[ match[i] ]]))
  pred.glm <- predict(m.glm, newdata = data.frame(x = datlist[[i]]$x))
  
  # random matching
  rid <- sample(1:length(modlist), 1)
  m.rand <- readRDS(paste0(path, "models_stan_eros/", modlist[[ rid ]]))
  pred.rand <- predict(m.rand, newdata = data.frame(x = datlist[[i]]$x))
  
  return(list(glm = pred.glm, rand = pred.rand))
}
stopCluster(cl)

errvec.glm.out <- rep(NA, length(modlist))
errvec.rand.out <- rep(NA, length(modlist))
for(i in 1:length(modlist)) {
  errvec.glm.out[i] <- mean(abs(datlist[[i]]$y - lout[[i]][[1]][,1]))
  errvec.rand.out[i] <- mean(abs(datlist[[i]]$y - lout[[i]][[2]][,1]))
}

plot(errvec.glm.out ~ mahdist, col = "gold", pch = 19, ylim = c(0, 12))
points(errvec.rand.out ~ mahdist, col = "magenta", pch = 19)

par(mfrow = c(2, 2))
boxplot(errvec.glm.out[ ,1:20], main = "GLM", ylab = "MAE (% cover)", xlab = "Time since fire (years)")

boxplot(errvec.glm.out[ ,1:20], main = "GLM", ylab = "Pr[ y_hat > y ]", xlab = "Time since fire (years)")
abline(h = .5, col = "red", lty = "dashed")}
# === end of draft code
#######################

#####################################################
library(leaflet)
library(viridisLite)
# --- plot fire polygons with leaflet ####
poly <- spTransform(tfires, CRS("+proj=longlat +datum=WGS84"))

errs <- errvec.glm.out[!is.na(errvec.glm.out)]

pal <- colorNumeric(viridis(50), domain = range(errs),
                    na.color = "transparent")
leaflet(poly[!is.na(errvec.glm.out), ]) %>%
  addPolygons(stroke = T, fillColor = pal, fillOpacity = 1) %>%
  addTiles(group = "OSM",
           options = providerTileOptions(minZoom = .1, maxZoom = 100)) %>%
  addProviderTiles("Esri.WorldTopoMap",    
                   group = "Topo") %>%
  addScaleBar(position = "bottomright")


leaflet() %>% 
  # addPolygons(data = poly, color = "black") %>%
  addPolygons(data = fires, color = "red") %>%
  addPolygons(data = fires, color = "blue") %>%
  addTiles(group = "OSM",
           options = providerTileOptions(minZoom = .1, maxZoom = 100)) %>%
  addProviderTiles("Esri.WorldTopoMap",    
                   group = "Topo") %>% 
  addScaleBar(position = "bottomright")


#### Export figures
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


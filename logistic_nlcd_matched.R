# === libraries + paths
pkgs <- c("brms", "matrixcalc", "tidyverse", "nlme", "lme4")
sapply(pkgs, require, character.only = T)

path <- "~/Desktop/SSM_PDE/nlcd_subset/"
years <- c(1985:2018)

# prepare intput data for non-spatial stan model
fires <- readOGR(paste0(path,"GB_wildfire_subset_1951_2018_ct_1_Rx0_blm_epsg4326.shp"))[1:10,]

sage <- readRDS(paste0(path,"sagelist_subs.rds"))
sagecrs <- sage[[1]]@crs

df <- spTransform(fires, sagecrs)
plot(df[c(4),])
plot(sage[[df$idvar[1]]][[1]],add=T)

#dfsubset = df[df$ha_clipped < 50 & df$Fire_Year > 1990,]
#sage1 = subset(sage, df$idvar %in% dfsubset$idvar)

mat <- as.data.frame(sage[[1]], xy = FALSE, na.rm = TRUE)
matplot(t(mat), type = "l", col = viridis(10), ylab = "Cover, %", xlab = "Time, years")
abline(v = df$FireYear[i] - 1984)
################## find the fire year based on the trajectories; Create offset variable
# n <- nrow(dfspat)
# offset <- rep(NA, n)
# yrs <- dfspat$Fire_Year - 1984
# sagemeans <- matrix(NA, nrow = n, ncol = 33)
# for(i in 1:n){
#   sagemeans[i, ] = cellStats(sage[[i]], mean)
#   f = which.min(diff(sagemeans[i, ])) + 1
#  if (yrs[i] == f - 1) {
#     offset[i] = 1
#   }  else if (yrs[i] == f - 2) {
#     offset[i] = 2
#   } else if (yrs[i] == f + 1) {
#     offset[i] = -1
#   } else {
#     offset[i] = 0
#   }
# }
# offset[250] <- -1
# offset[25] <- 1
# ind <- which(offset == 0)
# ind <- 1:nrow(dfspat)
# # plot the ts of average shrub values with proposed fire time point
# k = 9
# fig(index = dfspat$closeMatch)


# ====================================
# see pxlmatching.R 
tfires <- readRDS(paste0("~/Desktop/SSM_PDE/nlcd_subset/tfires.rds"))
tsage <- readRDS(paste0("~/Desktop/SSM_PDE/nlcd_subset/tsage.rds"))
tpxlcov <- readRDS(paste0("~/Desktop/SSM_PDE/nlcd_subset/tpxlcov.rds"))

kvec <- sapply(tpxlcov, function(x){ mean(x[['prefire']])})

# --- parametrize k prior
library(fitdistrplus)
klist <- list()
for(i in 1:length(tsage)) { klist[[i]] = tpxlcov[[i]]$prefire }
kprior <- unlist(klist)
fgamma <- fitdist( kprior[kprior > 0], "gamma" )
plot(density(kprior[kprior > 0]), lwd = 2)
curve(dgamma(x, shape = fgamma$estimate[1], rate = fgamma$estimate[2]), add = T, col ="brown", lwd = 3)

# === glmer 
mlist.glm <- list()
mlist0.glm <- list()
coefs <- list()
datlist.glm <- list()
tvec <- rep(NA, length(kvec))

for(i in 1:length(tsage)){
  print(i)
  
  # estimate growth and K parameters
  mat <- data.matrix(tsage[[i]][,c(tfires$FireYer[i]-1984):33])
  xpr <- ifelse(vec(mat[,-ncol(mat)]) == 0, .2, vec(mat[,-ncol(mat)]))
  ypr <- ifelse(vec(mat[,-1]) == 0, .2, vec(mat[,-1]))
  dat <- data.frame(y = log(ypr/xpr),
                    x = log(xpr),
                    cl = as.factor(rep(tpxlcov[[i]]$cluster, times = ncol(mat)-1)))
  temp <- lmer(y ~ (1|cl) + (0+x|cl), REML = FALSE, data = dat) #
  
  # estimate initial population size based the first 5 years of data post-fire
  mat <- tsage[[i]][,c(tfires$FireYer[i]-1984):33]
  T = ncol(mat)
  colnames(mat) <- 1:T-1
  mat$cl = as.factor(tpxlcov[[i]]$cluster)
  dat <- pivot_longer(mat, cols = c(1:T)) %>% 
    rename(t = name, cover = value) %>% 
    mutate(t = as.numeric(t), 
           cover = ifelse(cover == 0, runif(1e6, 0,1), cover)) %>%
    filter(t <= 5) %>% 
    select(cl, t, cover)
  temp_n0 <- glmer(cover ~ (1|cl) + (0+t|cl), family = Gamma(link = "log"), data = dat) # 

  # store models
  mlist.glm[[i]] = temp
  mlist0.glm[[i]] = temp_n0
  tvec[i] <- T
  datlist.glm[[i]] <- mat
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

# === in-sample predictions 
pred.glm <- list()
errvec.glm <- rep(NA, times = length(tsage))

for(j in 1:length(tsage)) {
  a <- coefs[[j]]$a
  b <- coefs[[j]]$b
  n0 <- coefs[[j]]$n0
  t <- 1:tvec[j]-1
  matpred <- gomp(a, b, n0, t)
  
  pred.glm[[j]] <- matpred
  a <- datlist.glm[[j]] %>% #select(-cl) #%>%
    group_by(cl) %>% 
    summarise_all(mean)
  errvec.glm[j] <- mean( abs(matpred - as.matrix(a[,-1])) )
}

errvec.glm[is.infinite(errvec.glm)] <- NA
plot(errvec.glm~kvec, xlab = "Pre-fire cover", ylab = "MAE", pch = 19)

# === out-of-sample errors; match sites at the fire level and pixel-to-cluster 
tdfEnv <- readRDS(paste0("~/Desktop/SSM_PDE/nlcd_subset/tdfEnv_covars.rds"))
tdfEnv$sm03mean <- ifelse(is.na(tdfEnv$sm03mean), mean(tdfEnv$sm03mean, na.rm = T), tdfEnv$sm03mean)
tdfEnv$sm03sd <-ifelse(is.na(tdfEnv$sm03sd), mean(tdfEnv$sm03sd, na.rm = T), tdfEnv$sm03sd)
tdfEnv$smAnnAnomsd <-ifelse(is.na(tdfEnv$smAnnAnomsd), mean(tdfEnv$smAnnAnomsd, na.rm = T), tdfEnv$smAnnAnomsd)

tdfEnv$u0 <- sapply(coefs, function(x) {mean(x$n0)} )

# --- site-to-site match
mah <- StatMatch::mahalanobis.dist(tdfEnv[,])
diag(mah) <- NA
match <- apply(mah, 1, which.min)
mahdist <- apply(mah, 1, min, na.rm = T)

# --- pixel-to-cluster match
pxl.dist <- function(dat, dat.m) {
  # finds a match for each focal pixel to a cluster in the reference site
  # - replaces cluster value for matched cluster
  # dat: focal pixel-level covariates
  # dat.m: reference pixel-level covariates
  df <- dat.m %>%
    group_by(cluster) %>%
    summarize_all(mean) %>% select(-cluster)
  mah.p <- StatMatch::mahalanobis.dist(dat[,-6], df)
  
  dat$cluster <- apply(mah.p, 1, which.min)
  return(dat)
}

errvec.glm.out <- rep(NA, length(tsage))
errvec.glm.out.rand <- rep(NA, length(tsage))

for(i in 1:length(tsage)){
  # glmer out-of-sample Mahalanobis nearest match
  j <- match[i]
  a <- coefs[[j]]$a
  b <- coefs[[j]]$b
  n0 <- coefs[[j]]$n0
  t <- 1:tvec[i]-1
  
  matpred <- as.data.frame(gomp(a, b, n0, t)) %>% mutate(cl = row_number())
  pred.m <- pxl.dist(tpxlcov[[i]], tpxlcov[[j]]) %>%
    select(cluster) %>% 
    left_join(matpred, by = c("cluster" = "cl")) %>%
    select(-cluster)
  
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
    select(cluster) %>% 
    left_join(matpred.r, by = c("cluster" = "cl")) %>%
    select(-cluster)
  
  obs <- datlist.glm[[i]] %>% select(-cl)
  
  errvec.glm.out[i] <- mean( as.matrix( abs(pred.m - obs)) )
  errvec.glm.out.rand[i] <- mean( as.matrix( abs(pred.r - obs)) )
}
errvec.glm.out[!is.finite(errvec.glm.out)] <- NA
errvec.glm.out.rand[!is.finite(errvec.glm.out.rand)] <- NA

plot(errvec.glm.out ~ mahdist, col = "purple", pch = 19, ylim = c(0, 12))
points(errvec.glm.out.rand ~ mahdist, col = "gold", pch = 19, ylim = c(0, 12))


plot(density(errvec.glm.out, na.rm = T), col = rgb(1,.5,1,1), xlab = "MAE")
lines(density(errvec.glm.out.rand, na.rm = T),col = rgb(0,.5,.1,1))

cat("In-sample MAE (%): ", mean(errvec.glm))
cat("Out-of-sample MAE (%): ", mean(errvec.glm.out))
cat("NULL out-of-sample MAE (%): ", mean(errvec.glm.out.rand))

##################################
# --- multiple matches averaged
mah <- StatMatch::mahalanobis.dist(tdfEnv[,])
diag(mah) <- NA
match <- apply(mah, 1, function(x) { which(x < min(x, na.rm = T)+1) } )

errvec.glm.out <- rep(NA, length(tsage))
errvec.glm.out.rand <- rep(NA, length(tsage))

for(i in 1:length(tsage)){
  # glmer out-of-sample Mahalanobis nearest match
  j <- match[[i]]
  amean <- mean( sapply(coefs[j], function(x) { x$a }) )
  bmean <- mean( sapply(coefs[j], function(x) { x$b }) )
  n0mean <- mean( sapply(coefs[j], function(x) { x$n0 }) )
  matpred <- matrix(NA, 1, tvec[i])
  for(k in 1:tvec[i]) { 
    matpred[,k] <- gomp(amean, bmean, n0mean, k-1)
  }
  # glmer out-of-sample, random match
  l <- sample(1:length(tsage), 1, replace = T)
  amean <- mean( sapply(coefs[l], function(x) { x$a }) )
  bmean <- mean( sapply(coefs[l], function(x) { x$b }) )
  n0mean <- mean( sapply(coefs[l], function(x) { x$n0 }) )
  matpred.rand <- matrix(NA, 1, tvec[i])
  for(k in 1:tvec[i]) { 
    matpred.rand[,k] <- gomp(amean, bmean, n0mean, k-1)
  }
  
  a <- datlist.glm[[i]] %>% select(-cl) %>%
    # group_by(cl) %>% 
    summarise_all(mean) 
  
  errvec.glm.out[i] <- mean(abs(matpred - as.matrix(a[,])))
  errvec.glm.out.rand[i] <- mean(abs(matpred.rand - as.matrix(a[,])))
}
errvec.glm.out[!is.finite(errvec.glm.out)] <- NA
errvec.glm.out.rand[!is.finite(errvec.glm.out.rand)] <- NA

plot(errvec.glm.out ~ mahdist, col = "purple", pch = 19, ylim = c(0, 12))
points(errvec.glm.out.rand ~ mahdist, col = "gold", pch = 19, ylim = c(0, 12))


plot(density(errvec.glm.out, na.rm = T), col = rgb(1,.5,1,1), xlab = "MAE")
lines(density(errvec.glm.out.rand, na.rm = T),col = rgb(0,.5,.1,1))

cat("In-sample MAE (%): ", mean(errvec.glm))
cat("Out-of-sample MAE (%): ", mean(errvec.glm.out))
cat("NULL out-of-sample MAE (%): ", mean(errvec.glm.out.rand))

################################################################################ 
res1 <- readRDS("~/Downloads/res1.rds")
ggplot(res1 , aes(y = value, x = as.factor(idvar))) + # geom_pointrange(aes(ymin = lower, ymax = upper), size = .1) +
  geom_boxplot() + 
  # scale_colour_viridis("Pre-disturbance abundance") +
  labs(y = "MAE", x = "ID") +  # geom_abline(intercept = 0, slope = 1) +
  theme_bw() + theme(axis.text.x = element_blank())

res1a <- res1[which(res1$mean < 5 & res1$k_mean > 10), ]

#####################################################
# --- plot fire polygons with leaflet ####
poly <- spTransform(df, CRS("+proj=longlat +datum=WGS84"))

# r1=sage[[ids1[1]]][[1]]
# r2=sage[[ids1[2]]][[1]]
#crs(r1) = sp::CRS("+init=epsg:3857")

pal <- colorNumeric(viridis(5), values(r1),
                    na.color = "transparent")
leaflet(poly[poly$Fire_Name == "Soda",]) %>%
  addPolygons(stroke = T, fillColor = "red", fillOpacity = 1) %>%
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


# === Export figures
# --- export leaflet/html visuals 
# library(htmlwidgets)
# saveWidget(m, file = paste0(pathfig, "figmap.html"))
# library(mapview)
# mapshot(m, file = paste0(pathfig, "figmap1.pdf"))

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


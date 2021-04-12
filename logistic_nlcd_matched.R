# === libraries + paths
pkgs <- c("brms", "matrixcalc", "reshape2", "nlme", "lme4")
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
plot(density(kprior), lwd = 2)
curve(dgamma(x, shape = fgamma$estimate[1], rate = fgamma$estimate[2]), add = T, col ="brown", lwd = 3)

# --- Poisson GLM on the log-log scale
# temp <- brm(y ~ 1 + x  + offset(x), family = "poisson", data=dat, chains = 1, iter = 100)

# === glmer predictions
errvec.glm <- rep(NA, times = length(tsage))
mlist.glm <- list()
datlist.glm <- list()
for(i in 1:length(tsage)){
  print(i)
  mat <- data.matrix(tsage[[i]][,c(tfires$FireYer[i]-1984):33])
  xpr <- vec(mat[,-ncol(mat)])
  dat <- data.frame(y = vec(mat[,-1]), 
                    x = log(ifelse(xpr == 0, .25, xpr)),
                    cl = as.factor(rep(tpxlcov[[i]]$cluster, times = ncol(mat)-1)))
  #temp <- glmer(y ~ (1|cl) + x + (x|cl) + offset(x), family = "poisson", data = dat)
  temp <- brm(y ~ (1|cl) + (x|cl)  + offset(x), family = "poisson", data=dat, chains = 1, iter = 100)
  pred <- exp(predict(temp))
  errvec.glm[i] <- mean(abs(pred - dat$y))
  
  mlist.glm[[i]] = temp
  datlist.glm[[i]] = dat
}
plot(errvec.glm ~ kvec, pch = 19)

# exact Logistic solution
# K * N0 * exp(r * t) / (K + N0 * (exp(r * t) - 1))

mat <- tsage[[i]][,c(tfires$FireYer[i]-1984):33]
colnames(mat) <- 1:ncol(mat)
mat$n0 = as.numeric(ifelse(mat[,1] == 0,  0.01, mat[,1]))
mat$cl = as.factor(tpxlcov[[i]]$cluster)
dat <- melt(mat, id.vars = c("n0", "cl"))
names(dat) <- c("n0", "cl", "t", "cover")
dat$t <- as.numeric(dat$t)

lgnl <- bf(cover ~ u0*k/( u0 + (k - u0)*exp(-r*t) ),
             r ~ 1|cl, k ~ 1|cl, u0 ~ 1, nl = TRUE)
nlprior <- c(set_prior("gamma(4.085748, 0.2505014)", lb = 0, nlpar = "k"),
             set_prior("gamma(3, 3)", lb = 0, nlpar = "r"),
             set_prior("exponential(1)", lb = 0, nlpar = "u0"))
fit1 <- brm(formula = lgnl, data = dat, family = gaussian(),
            prior = nlprior, 
            control = list(adapt_delta = 0.9),
            chains = 1, iter = 200)
stanfit <- fit1$fit

fit2 <- brm(formula = lgnl, data = dat, family = gaussian(),
            prior = nlprior, 
            control = list(adapt_delta = 0.9),
            chains = 1, iter = 200)
plot( conditional_effects(fit1), points = T)


errvec.lg <- rep(NA, times = length(tsage))
mlist.lg = list()
datlist.lg = list()
for(i in 1:length(tsage)){
  print(i)
  mat <- tsage[[i]][,c(tfires$FireYer[i]-1984):33]
  colnames(mat) <- 1:ncol(mat)
  mat$n0 <- mat[, 1]
  mat$cl <- as.factor(tpxlcov[[i]]$cluster)
  dat <- melt(mat, id.vars = c("n0", "cl"))
  names(dat) <- c("n0", "cl", "t", "cover")
  dat$t <- as.numeric(dat$t)

  m1 <- nlme(cover ~ u0*k/( u0 + (k - u0)*exp(-r*t)),
             data = dat, 
             fixed = k + u0 + r ~ 1,
             random =  r + k ~ 1,
             groups = ~cl,
             start = c(k = 20, u0 = 1, r = 1))

  pred <- predict(m1)
  errvec.lg[i] <- mean(abs(pred - dat$cover))
  
  mlist.lg[[i]] <- m1
  datlist.lg[[i]] <- dat
}
plot(errvec.lg ~ kvec, pch = 19, col = rgb(.5,0,.5,1))

plot(density(errvec.glm), xlim = c(.5, 5), col = "purple", lwd = 2)
lines(density(errvec.lg), lwd = 2, col = "brown")

# --- in-sample error estimates
err.yr.glm <- matrix(NA, 12, 33)
err.yr.lg <- matrix(NA, 12, 33)
prob.yr.glm <- matrix(NA, 12, 33)
prob.yr.lg <- matrix(NA, 12, 33)
for(i in 1:length(tsage)) {
  pred1 = exp(predict(mlist.glm[[i]]))
  pred2 = predict(mlist.lg[[i]])
  dat1 = datlist.glm[[i]]
  dat2 = datlist.lg[[i]]
  names(pred1) = dat2$t[dat2$t > 1]
  dat1$t = dat2$t[dat2$t > 1]
  for(t in 2:33) {
    # mean absolute error
    err.yr.glm[i,t] = mean(abs(pred1[dat1$t == t] - dat1$y[dat1$t == t]))
    err.yr.lg[i,t-1] = mean(abs(pred2[dat2$t == t-1] - dat2$cover[dat2$t == t-1]))
    # probability of over-predicting
    # prob.yr.glm[i,t] = sum(pred1[dat1$t == t] > dat1$y[dat1$t == t]) / length(dat1$y[dat1$t == t])
    # prob.yr.lg[i,t-1] <- sum(pred2[dat2$t == t-1] > dat2$cover[dat2$t == t-1]) / length(dat2$cover[dat2$t == t-1])
    # average residuals by year
    prob.yr.glm[i,t] <- mean(pred1[dat1$t == t] - dat1$y[dat1$t == t]) 
    prob.yr.lg[i,t-1] <- mean(pred2[dat2$t == t-1] - dat2$cover[dat2$t == t-1])
  }
}

par(mfrow = c(2, 2))
boxplot(err.yr.glm[ ,1:20], main = "GLM", ylab = "MAE (% cover)", xlab = "Time since fire (years)")
boxplot(err.yr.lg[ ,1:20], main = "Logistic growth", ylab = "MAE (% cover)", xlab = "Time since fire (years)")

boxplot(prob.yr.glm[ ,1:20], main = "GLM", ylab = "Residuals [y_hat - y]", xlab = "Time since fire (years)")
abline(h = .5, col = "red", lty = "dashed")
boxplot(prob.yr.lg[ ,1:20], main = "Logistic growth", ylab = "Residuals [ y_hat - y ]", xlab = "Time since fire (years)")
abline(h = .5, col = "red", lty = "dashed")

# --- out-of-sample errors; match sites at the fire level
tdfEnv <- readRDS(paste0("~/Desktop/SSM_PDE/nlcd_subset/tdfEnv_covars.rds"))

tdfEnv$u0 <- sapply(mlist.lg, function(x) {x$coefficients$fixed['u0']} )
  
mah <- StatMatch::mahalanobis.dist(tdfEnv)
diag(mah) <- NA
match <- apply(mah, 1, which.min)
mahdist <- apply(mah, 1, min, na.rm = T)

errvec.glm.out <- rep(NA, length(tsage))
errvec.lg.out <- rep(NA, length(tsage))
errvec.rand.out <- rep(NA, length(tsage))
for(i in 1:length(tsage)){
  # glmer out-of-sample predictions
  pred1 <- exp(predict(mlist.glm[[ match[i] ]], newdata = data.frame(x = datlist.glm[[i]]$x, cl = datlist.glm[[i]]$cl)))
  errvec.glm.out[i] <- mean(abs(datlist.glm[[i]]$y - pred1))
  # # non-linear out-of-sample predictions
  # pred2 <- predict(mlist.lg[[ match[i] ]], newdata = data.frame(t = datlist.lg[[i]]$t, cl = datlist.lg[[i]]$cl))
  # errvec.lg.out[i] <- mean(abs(datlist.lg[[i]]$cover - pred2))
  # random matching
  rid <- sample(1:12, 1)
  pred3 <- predict(mlist.lg[[rid]], newdata = data.frame(t = datlist.lg[[i]]$t, cl = datlist.lg[[i]]$cl))
  errvec.rand.out[i] <- mean(abs(datlist.lg[[i]]$cover - pred3))
}

plot(errvec.glm.out ~ mahdist, col = "purple", pch = 19, ylim = c(0, 12))
points(errvec.lg.out ~ mahdist, col = "green", pch = 19)
points(errvec.rand.out ~ mahdist, col = "magenta", pch = 19)

par(mfrow = c(2, 2))
boxplot(err.yr.glm[ ,1:20], main = "GLM", ylab = "MAE (% cover)", xlab = "Time since fire (years)")
boxplot(err.yr.lg[ ,1:20], main = "Logistic growth", ylab = "MAE (% cover)", xlab = "Time since fire (years)")

boxplot(prob.yr.glm[ ,1:20], main = "GLM", ylab = "Pr[ y_hat > y ]", xlab = "Time since fire (years)")
abline(h = .5, col = "red", lty = "dashed")
boxplot(prob.yr.lg[ ,1:20], main = "Logistic growth", ylab = "Pr[ y_hat > y ]", xlab = "Time since fire (years)")
abline(h = .5, col = "red", lty = "dashed")





################################## 
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


#### Export figures
# ----
library(export)
pathfig <- "~/Downloads/"
graph2eps(file = paste0(pathfig, "figh2.eps"), height = 6, width = 6)

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


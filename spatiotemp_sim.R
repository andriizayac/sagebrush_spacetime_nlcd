 ###### Spatio temporal data simulation
library(raster)
library(rstan)
library(ggplot2)
library(dplyr)
library(viridisLite)
library(reshape2)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

neighborhood <- function(raster){
  nn <- matrix(NA,length(raster[]),4)
  for(i in 1:dim(nn)[1]){
    loc <- adjacent(raster,i)[,2]
    ln <- loc[which((loc+1)==i)]
    rn <- loc[which((loc-1)==i)]
    bn <- loc[which((loc-dim(raster)[2])==i)]
    tn <- loc[which((loc+dim(raster)[2])==i)]
    nn[i,1] <- ifelse(length(ln)>0,ln,0)
    nn[i,2] <- ifelse(length(rn)>0,rn,0)
    nn[i,3] <- ifelse(length(tn)>0,tn,0)
    nn[i,4] <- ifelse(length(bn)>0,bn,0)
  }
  nn
}
neighborhoodMat <- function(raster){
  N <- length(raster[])
  nn <- matrix(0, N, N)
  diag(nn) <- 1:N
  for(i in 1:dim(nn)[1]){
    loc <- adjacent(raster, i)[, 2]
    ln <- loc[which((loc + 1) == i)]
    rn <- loc[which((loc - 1) == i)]
    bn <- loc[which((loc - dim(raster)[2]) == i)]
    tn <- loc[which((loc + dim(raster)[2]) == i)]
    if (length(ln) > 0) { nn[i, i+1] = ln } 
    if (length(rn) > 0) { nn[i, i-1] = rn }
    nn[i,3] <- ifelse(length(tn) > 0, tn, 0)
    nn[i,4] <- ifelse(length(bn) > 0, bn, 0)
  }
  nn
}
wrap <- function(x, poly){
  row = dim(poly)[1]
  col = dim(poly)[2]
  mat=matrix(x, row, col, byrow=TRUE)
  mat
}

### - a function of discretized Fickian diffusion
makeH <- function(NN, Dvec, dx, dy){
  H <- matrix(0, dim(NN)[1], dim(NN)[1])
  diag(H) <- Dvec * (1 / dx^2 + 1 / dy^2)
  for(i in 1:dim(H)[1]){
    H[i,i] = Dvec[i]*(1/dx^2 + 1/dy^2)
    if(sum(NN[i, ] != 0) == 4){
      if(NN[i, 1] > 0){H[i, NN[i, 1]] <- (Dvec[i] - (Dvec[NN[i, 2]] - Dvec[NN[i, 1]]) / 4) / (dx^2)}
      if(NN[i, 2] > 0){H[i, NN[i, 2]] <- (Dvec[i] + (Dvec[NN[i, 2]] - Dvec[NN[i, 1]]) / 4) / (dx^2)}
      if(NN[i, 3] > 0){H[i, NN[i, 3]] <- (Dvec[i] + (Dvec[NN[i, 4]] - Dvec[NN[i, 3]]) / 4) / (dy^2)}
      if(NN[i, 4] > 0){H[i, NN[i, 4]] <- (Dvec[i] - (Dvec[NN[i, 4]] - Dvec[NN[i, 3]]) / 4) / (dy^2)}
    }}
  H
}
makeG <- function(NN, uvec , dx, dy){
  G = matrix(0,dim(NN)[1],dim(NN)[1])
  for(i in 1:dim(H)[1]){
    if(sum(NN[i,]!=0) == 4){
      G[i,i] <- (uvec[NN[i,2]] - 2*uvec[i] + uvec[NN[i,1]])/(2*dx^2) + (uvec[NN[i,4]] - 2*uvec[i] + uvec[NN[i,3]])/(2*dy^2)
      if(NN[i,1]>0){H[i,NN[i,1]] <- -(uvec[NN[i,2]] - uvec[NN[i,1]])/(4*dx^2)}
      if(NN[i,2]>0){H[i,NN[i,2]] <- (uvec[NN[i,2]] - uvec[NN[i,1]])/(4*dx^2)}
      if(NN[i,3]>0){H[i,NN[i,3]] <- -(uvec[NN[i,3]] - uvec[NN[i,3]])/(4*dy^2)}
      if(NN[i,4]>0){H[i,NN[i,4]] <- (uvec[NN[i,4]] - uvec[NN[i,3]])/(4*dy^2)}}
  }
  G
}

### - a function of discretized Plain diffusion with varying D
makeGplain <- function(NN, uvec , dx, dy, dt){
  G <- matrix(0, dim(NN)[1], dim(NN)[1])
  for(i in 1:dim(G)[1]){
    G[i,i] <- -4 * uvec[i]
    if(sum(NN[i, ] != 0) == 4){
      if(NN[i, 1] > 0){G[i, NN[i, 1]] <- uvec[NN[i,1]]}
      if(NN[i, 2] > 0){G[i, NN[i, 2]] <- uvec[NN[i,2]]}
      if(NN[i, 3] > 0){G[i, NN[i, 3]] <- uvec[NN[i,3]]}
      if(NN[i, 4] > 0){G[i, NN[i, 4]] <- uvec[NN[i,4]]}}
  }
  G *(dt / (dx^2 + dy^2))
}
### - a neighborhood matrix 0/1
#makeNNBmatrix <- function()

# --- load data and transform to N x T formal
# years <- c(1985:2017)
# path <- "~/MEGA/landsat_backintime/"
# file <- "BogusCreek_1986_clipped_30_30.tif" #1986
# ras <- stack(paste0(path,file))
# poly <- subset(ras,1:14)
# values(poly)[values(poly) == 101] <- NA
# df1 <- as.data.frame(poly, xy = TRUE, na.rm = TRUE)
# names(df1) <- c('x', 'y', paste0('year', years[1:14]))

# === a setup for deterministic simulation
n <- 30
T <-22
matsim <- matrix(runif(n^2, 0.001, .002), n, n)
NN0 <- neighborhood(raster(matsim))
bndIdx = which(apply(NN0, 1, function(x) {sum(x == 0) > 0 }))
NN <- NN0
# make diffusion stay within the spatial field
for(i in 1:nrow(NN)){
  for(j in 1:4){
    if(NN[i, j] == 0) {
      NN[i, j] <- i
    }
  }
}
adjMatrix = landscapemetrics::get_adjacencies(raster(matsim), neighbourhood = 4, what = "full")

# mesh parameters
dx <- 30
dy <- 30
dt <- 1

# demographic parameters
r <- .7
k <- 1
D <- rep(100, n^2)

# initiate the simulation
# var <- sample(n, n/3)
matsim[var, var] <- k/2 #rnorm(25, 2, .1)
# matsim[1:2, ] <- k/2 #rnorm(25, 2, .1)
ysim <- matrix(0, n^2, T) #runif(n^2, 0, k/3)
ysim[, 1] <- as.numeric(matsim)
ysim <- ysim

ysim <- simsage(ysim = ysim, r = r, k = k, D = D, 
               dx = dx, dy = dy, dt = dt, NN = NN,
               spat = FALSE)

matplot(t(ysim[,]) , type = "l" , ylab = "Cover, %", xlab = "Time, years", col = viridis(n))
# plot raster time-series
par(mfrow = c(3, 3), mar = c(1, 2.5, 1, 2))
for(i in 1:T){
  ymat = wrap(ysim[, i], # t(mat[i, ]), #
              raster(matsim))
  plot(raster(ymat), col = viridis(25))
  rastlist[[i]] <- raster(ymat)
}
dev.off()
raststack <- stack(rastlist)
animate(raststack, pause=1, main = c(1:9), zlim = c(0,1), maxpixels=900, n=1, col = viridis(10))
saveGIF(for(i in 1:T){raster::plot(raststack[[i]], col = viridis(20), main = paste("Time:", i))}, 
        movie.name = "~/Downloads/animation2.gif", interval=1)


# === stability analysis ####
dt <- rep(1, 9) #c(.1, .5, .75, 1, 10, 100, 150, 200, 365) #seq(.6, 1.4, l = 9)
Dvec <-  seq(1, 500, l = 9) # rep(1, 9) #
D <- matrix(NA, n^2, T)
par(mfrow = c(3, 3), mar = c(2, 1, 2.5, 1))
for(i in 1:length(Dvec)){
  D[, i] = Dvec[i]
  ysim <- simsage(ysim = ysim, r = .9, k = k, D = D[, i], 
                  dx = dx, dy = dy, dt = dt[i], NN = NN,
                  spat = TRUE)
  CFL = round((mean(Dvec[i]) * dt[i]) / dx, 4)
  matplot(t(ysim[-bndIdx,]), type="l" , ylab="Cover, %", xlab="Time, years", col = viridis(n), 
          main = paste("dt =", dt[i], ", CFL =", CFL))
}

# === a setup for stochastic state-space simulation
n <- 15
T <- 12
matsim <- matrix(runif(n^2, 2, 4), n, n)
NN0 <- neighborhood(raster(matsim))
bndIdx = which(apply(NN0, 1, function(x) {sum(x == 0) > 0 }))
NN <- NN0
# make diffusion stay within the spatial field
for(i in 1:nrow(NN)){
  for(j in 1:4){
    if(NN[i, j] == 0) {
      NN[i, j] <- i
    }
  }
}

# mesh parameters
dx <- 30
dy <- 30
dt <- 1

# demographic parameters
r <- .7
D <- rep(100, n^2) # rnorm(n^2, 100, .5)

# errors
gamma <- .01
eta <- seq(.03, .05, l = T)

# initiate the simulation
zinit <- matrix(runif(n^2, .1, .11), n, n)
zmat <- matrix(0, n^2, T+1)
zmat[, 1] <- as.numeric(zinit)
ymat <- matrix(0, n^2, T) 

ysim <- simsage_statespace(ysim = zmat, ymat = ymat, r = r, k = k, D = D, 
                dx = dx, dy = dy, dt = dt, NN = NN, spat = FALSE,
                gamma, eta)

matplot(t(ysim[[2]]), type = "l" , ylab = "Cover, %", xlab = "Time, years", col = viridis(n))
par(mfrow = c(3, 3), mar = c(1, 1, 1, 1))
for(i in 1:T){
  ymat = wrap(t(ysim[[1]][, i]), #ysim[, i], # t(mat[i, ]), #
              raster(matsim))
  plot(raster(ymat))
  #plot(poly[[i]])
}
dev.off()


# === linearized models
dat <- melt(ysim)
names(dat)=c("id","year","cover")
plot(1,type="n",xlim=c(0,max(dat$year)),ylim=c(-2,4))
for(i in 1:nrow(ysim)){
  points(lcover~year,data=dat[dat$id==i,],type="l")
}
dat$lcover = log(dat$cover)
lm(dat$lcover~log(dat$year))
curve(-.75+1.5*log(x),0,12,col="red",lwd=5,add=TRUE)
points(1,log(mean(ysim[,1])),col="green",cex=2,pch=19)

# === stan model 
dlist <- list(N = n^2, T = T, 
              y_mat = t(ysim[[2]]), y0_tr = as.numeric(zinit),#ysim[[2]][, 1], 
              NN = NN)
library(rstan)
m <- stan(file = "rdiff_discrete_plain_statespace.stan",
              data = dlist, chains = 3, iter = 300, warmup = 275)
print(m, pars = c("r"))
plot(m, pars = c("eta"))
# plot(m, pars = c("Delta"))

mat <- with(rstan::extract(m), Z)
mat <- apply(mat, c(2,3), mean)
init <- with(extract(m), z_init)
init <- apply(init, 2, mean)

matplot(mat[-1, ], type = "l" , ylab = "Cover, %", xlab = "Time, years", col = viridis(n))
matplot(t(ysim[[1]]) * k, type = "l" , ylab = "Cover, %", xlab = "Time, years", col = viridis(n))
plot(mat[-1, ] ~ t(ysim[[1]]))
plot(init ~ as.numeric(zinit))
abline(0, 1)

vec <- matrix(1:n^2, n, n)
vec <- ifelse(vec %in% bndIdx, 4, 5)
plot((mat[-1, ] *khat) ~ t(ysim[[2]]), pch = 19)
abline(0, 1, col = "gold", lty = "dashed", lwd = 3)


par(mfrow = c(3, 3), mar = c(1, 1, 1, 1))
for(i in 1:T){
  mat1 <- t(mat[-1, ])
  a <- abs(mat1[i, ] - t(ysim[[1]])[i, ])
  plot(raster(wrap(a, raster(matsim))), col = viridis(10))
}

#### posterior predictions
genPred <- function(model) {
  mpost <- rstan::extract(model, data)
  r <- with(mpost, r)
  k_mu <- with(mpost, k_mu)
  z_init <- with(mpost, z_init)
  tj <- array(NA, dim = c(length(r), data$T,data$N))
  for(i in 1:length(r)) {
    for(j in 1:data$T){
    }
  }
}



adjacencies <- raster::adjacent(raster(matsim), 1:raster::ncell(raster(matsim)), 4, pairs=TRUE)
table(raster(matsim)[adjacencies[,1]], raster(matsim)[adjacencies[,2]])

a = matrix(1:9, 3, 3)
arast = raster(a)
c = adjacent(arast, 1:ncell(arast), 4, pairs = T, include = TRUE)
adj = table(arast[c[, 1]], arast[c[, 2]])


### === Gompertz model


#### Export figures ####
library(export)
pathfig <- "~/Downloads/"
graph2eps(file = paste0(pathfig, "fig_3.eps"), height = 3, width = 3.5)

 ###### Spatio temporal data simulation
library(raster)
library(ggplot2)
library(dplyr)
library(viridisLite)
library(reshape2)

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
wrap = function(x,poly){
  row=dim(poly)[1]
  col=dim(poly)[2]
  mat=matrix(x, row,col, byrow=TRUE)
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

# --- load data and transform to N x T formal
# years <- c(1985:2017)
# path <- "~/MEGA/landsat_backintime/"
# file <- "BogusCreek_1986_clipped_30_30.tif" #1986
# ras <- stack(paste0(path,file))
# poly <- subset(ras,1:14)
# values(poly)[values(poly) == 101] <- NA
# df1 <- as.data.frame(poly, xy = TRUE, na.rm = TRUE)
# names(df1) <- c('x', 'y', paste0('year', years[1:14]))

n <- 20
T <- 10
matsim <- matrix(0, n, n)
NN <- neighborhood(raster(matsim))
bndIdx = which(apply(NN, 1, function(x) {sum(x == 0) == 1 | sum(x == 0) == 1}))

ymat <- as.data.frame(matsim, xy = FALSE, na.rm = TRUE)

dx <- 30
dy <- 30
dt <- 1

r <- .4
k <- 15
D <- rep(50, n^2)
G <- makeGplain(NN, y0 , dx, dy, dt)

# initiate the simulation
var <- sample(n, n/3)
matsim[var, var] <- k
ysim <- matrix(0.1, n^2, T)
ysim[, 1] <- as.numeric(matsim)
ysim <- ysim/k

simsage <- function(ysim, r, k, D, dx, dy, dt, NN, spat = TRUE) {
  if(sum((ysim[, 1] > k)) > 0){
    return("error")
  } else{
    for(t in 2:ncol(ysim)){
      print(t)
      if(spat == TRUE){
        G = makeGplain(NN, ysim[, t-1], dx, dy, dt)
        ysim[, t] =  ysim[, t-1] + G %*% D + r * ysim[, t-1] * (1 -  ysim[, t-1]) 
      } else {
        ysim[, t] =  ysim[, t-1] + r * ysim[, t-1] * (1 -  ysim[, t-1])
      }
    }
  return(ysim)
  }
}

ysim <- simsage(ysim = ysim, r = r, k = k, D = D, 
               dx = dx, dy = dy, dt = dt, NN = NN,
               spat = TRUE)

matplot(t(ysim)*k, type="l" , ylab="Cover, %", xlab="Time, years", col = viridis(n))
# --- plot raster time-series
par(mfrow=c(3,3),mar=c(1,1,1,1))
for(i in 1:T){
  ymat = wrap(ysim[, i], 
              raster(matsim))
  plot(raster(ymat))
  #plot(poly[[i]])
}
dev.off()
# ---------------------------------------------------
Dvec <- seq(0, 300, l = 9)
D <- matrix(NA, n^2, T)
par(mfrow = c(3, 3), mar = c(1, 1, 1, 1))
for(i in 1:length(Dvec)){
  D[, i] = Dvec[i]
  ysim <- simsage(ysim = ysim, r = .1, k = k, D = D[, i], 
                  dx = dx, dy = dy, dt = dt, NN = NN,
                  spat = TRUE)
  matplot(t(ysim[-bndIdx,]), type="l" , ylab="Cover, %", xlab="Time, years", col = viridis(n))
}

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


dlist <- list(N = n^2, T = T, N_pred = n^2, T_pred = T,
              y_mat = t(ysim)*k, y0_tr = ysim[, 1], y0_pred = ysim[, 1],
              y_k_prior = rep(k, n^2), y_k_pred = rep(k, n^2),
              NN = NN)
library(rstan)
m <- stan(file = "rdiff_discrete_plain_statespace.stan",
              data = dlist, chains = 2, iter = 100, warmup = 90)
plot(m, pars = c("D"))

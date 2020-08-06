library(rstan)
library(raster)
library(rgdal)
library(leaflet)
library(viridis)

# === allow multi-core stan sa,mpling
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# === calculate Mean Absoute Error between x and y
mae=function(x,y){sum(abs(x-mean(y)))/length(x)}
# plot pre- and post- fire rasters (based on the dfspat - matched data frame) 
plotfire=function(id=NULL,offyr=offset,focal=TRUE,pol=dfspat,rast=sage){
  par(mfrow=c(1,2))
  i = ifelse(focal==T, dfspat$idvar[id], dfspat$closeMatch[id])
  k = dfspat$Fire_Year[i] - 1984
  plot(sage[[i]][[k - 1 + offyr[i]]], main=paste0(k -1 +1984,": pre"))
  plot(sage[[i]][[k + offset[i]]], main=paste0(k + 1984,": post"))
}
# prepare intput data for non-spatial stan model
dataprep <- function(n=NULL,pol=NULL,rast=NULL,years=NULL,yroffset=NULL,
                    seed=123,trsubset=.5,predsubset=1,threshold=FALSE,thresholdval=NULL){
  # n: #row in the matched polygon data frame
  # pol: polygon dataframe including target and ref sites
  # rast: list of rasters corresponding to each target/ref polygon (each has 33 layers=years)
  # years: vector of years - 1985:2017
  # yroffset: mismatch b/w fire year and NLCD fire; df w/ 2 cols: target, ref
  # trsubset: proportion to subset from the full dataset as training sample
  # predsubset: proportion to predict over in the test dataset
  # threshold: value [0,1] in shrub cover to include as a training dataset.Values closer to one will allow data points that introduce bias into 'r' parameter 
  
  set.seed(seed)
  n=n
  ids=c(pol$idvar[n],pol$closeMatch[n])
  fyear = c(pol$Fire_Year[ids[1]],pol$Fire_Year[ids[2]])
  if(fyear[1]<fyear[2]){
    # subset post-disturbance trajectories
    dftr = as.data.frame(rast[[ids[1]]],xy=F,na.rm = T)[, years >= fyear[1] + yroffset[ids[1]] ]
    dftest = as.data.frame(rast[[ids[2]]],xy=F,na.rm = T)[, years >= fyear[2] + yroffset[ids[2]] ]
    # create training indices
    tr=sample(1:nrow(dftr), round(nrow(dftr) * trsubset))
    tr_pred=sample(1:nrow(dftest), round(nrow(dftest) * predsubset))
    # create training and test datasets, create pre-disutrbance means
    y0_train = dftr[tr,1]
    if(sum(years < fyear[1] - 1) == 1) {
      y_k_prior = as.numeric(as.data.frame(rast[[ids[1]]],xy=F,na.rm = T)[tr, 1:sum(years < fyear[1] - 1)])
    } else {
      y_k_prior = as.numeric(apply(as.data.frame(rast[[ids[1]]],xy=F,na.rm = T)[tr, 1:sum(years < fyear[1] - 1)], 1, mean))
    }
    y0_predict = dftest[tr_pred,1]
    y_k_pred = as.numeric(apply(as.data.frame(rast[[ids[2]]],xy=F,na.rm = T)[tr_pred, 1:sum(years < fyear[2] - 1)], 1, mean))
    # standardize initial values to [0,1] for non-dimentional model
    y0_tr = ifelse(y0_train < y_k_prior, y0_train / y_k_prior, y0_train / y0_train)
    y0_pred = ifelse(y0_predict < y_k_pred, y0_predict / y_k_pred, y0_predict / y0_predict)
    if(threshold==TRUE){
      y0_train = ifelse(is.finite(y0_tr),y0_tr,0)
      idx = which(y0_train < thresholdval)
      while(length(idx) < 100) {
        thresholdval = thresholdval + .1
        idx = which(y0_train < thresholdval)
      }
      y_matrix = as.matrix(dftr[tr, ])[idx, ]
      
      dlist=list(N = nrow(y_matrix),
                 T = ncol(dftr),
                 # training data
                 y_mat=t(y_matrix),
                 y_k_prior = y_k_prior[idx],
                 y0_tr = y0_train[idx],
                 # predictions data
                 y_k_pred = y_k_pred,
                 y0_pred = ifelse(is.finite(y0_pred),y0_pred,0),
                 N_pred = length(tr_pred),
                 T_pred = ncol(dftest),
                 # record threshold used
                 cover_threshold = thresholdval,
                 # input data sets
                 dftr=as.data.frame(y_matrix),
                 dftest=dftest[tr_pred,])
    } else {
    dlist=list(N = length(tr),
               T = ncol(dftr),
               # training data sets
               y_mat=t(as.matrix(dftr[tr,])),
               y_k_prior = y_k_prior,
               y0_tr = ifelse(is.finite(y0_tr),y0_tr,0),
               # predictions data
               y_k_pred = y_k_pred,
               y0_pred = ifelse(is.finite(y0_pred),y0_pred,0),
               N_pred = length(tr_pred),
               T_pred = ncol(dftest),
               # input data sets
               dftr=dftr[tr,],
               dftest=dftest[tr_pred,])
    }
  } else {
    # subset post-disturbance trajectories
    dftr = as.data.frame(rast[[ids[2]]],xy=F,na.rm = T)[, years>=fyear[2] + yroffset[ids[2]] ]
    dftest = as.data.frame(rast[[ids[1]]],xy=F,na.rm = T)[, years>=fyear[1] + yroffset[ids[1]] ]
    # create training indices
    tr=sample(1:nrow(dftr),round(nrow(dftr)*trsubset))
    tr_pred=sample(1:nrow(dftest),round(nrow(dftest)*predsubset))
    # create training and test datasets, create pre-disutrbance means
    y0_train = dftr[tr,1]
    if(sum(years < fyear[2] - 1) == 1){
      y_k_prior = as.numeric(as.data.frame(rast[[ids[2]]],xy=F,na.rm = T)[tr, 1:sum(years < fyear[2] - 1)])
    } else {
      y_k_prior = as.numeric(apply(as.data.frame(rast[[ids[2]]],xy=F,na.rm = T)[tr, 1:sum(years < fyear[2] - 1)], 1, mean))
    }
    y0_predict = dftest[tr_pred,1]
    y_k_pred = as.numeric(apply(as.data.frame(rast[[ids[1]]],xy=F,na.rm = T)[tr_pred, 1:sum(years < fyear[1] - 1)], 1, mean))
    # standardize initial values to [0,1] for non-dimentional model
    y0_tr = ifelse(y0_train < y_k_prior, y0_train / y_k_prior, y0_train / y0_train)
    y0_pred = ifelse(y0_predict < y_k_pred, y0_predict / y_k_pred, y0_predict / y0_predict)
    if(threshold==TRUE){
      y0_train = ifelse(is.finite(y0_tr), y0_tr, 0)
      idx=which(y0_train < thresholdval)
      while(length(idx) < 100) {
        thresholdval = thresholdval + .1
        idx=which(y0_train < thresholdval)
      }
      y_matrix = as.matrix(dftr[tr, ])[idx, ]
      
      dlist=list(N = nrow(y_matrix),
                 T = ncol(dftr),
                 # training data
                 y_mat=t(y_matrix),
                 y_k_prior = y_k_prior[idx],
                 y0_tr = y0_train[idx],
                 # predictions data
                 y_k_pred = y_k_pred,
                 y0_pred = ifelse(is.finite(y0_pred),y0_pred,0),
                 N_pred = length(tr_pred),
                 T_pred = ncol(dftest),
                 # record threshold used
                 cover_threshold = thresholdval,
                 # input data sets
                 dftr=as.data.frame(y_matrix),
                 dftest=dftest[tr_pred,])
    } else {
    dlist=list(N = length(tr),
               T = ncol(dftr),
               # training data
               y_mat=t(as.matrix(dftr[tr,])),
               y_k_prior = y_k_prior,
               y0_tr = ifelse(is.finite(y0_tr),y0_tr,0),
               # predictions data
               y_k_pred = y_k_pred,
               y0_pred = ifelse(is.finite(y0_pred),y0_pred,0),
               N_pred = length(tr_pred),
               T_pred = ncol(dftest),
               # input data sets
               dftr=dftr[tr,],
               dftest=dftest[tr_pred,])
    }
  }
  return(dlist)
}


years=c(1985:2017)
path = "C:/Users/CaughlinLab/Desktop/Landsat_eros/"

sage=readRDS(paste0(path,"usgs_sagebrush/sagelist_finite.rds"))
dfspat = readRDS(paste0(path,"mtch_polygons.rds"))

df=spTransform(dfspat,sage[[1]][[1]]@crs)
plot(df[c(2,43),])
plot(sage[[df$idvar[2]]][[1]],add=T)
plot(sage[[df$idvar[43]]][[1]],add=T)

#dfsubset = df[df$ha_clipped < 50 & df$Fire_Year > 1990,]
#sage1 = subset(sage, df$idvar %in% dfsubset$idvar)

plot(sage[[1]][[1]])
plot(dfspat[1,],add=T)
################## find the fire year based on the trajectories; Create offset variable
n = nrow(dfspat)
offset = rep(NA, n)
yrs = dfspat$Fire_Year - 1984
sagemeans = matrix(NA, nrow = n, ncol = 33)
for(i in 1:n){
  sagemeans[i, ] = cellStats(sage[[i]], mean)
  f = which.min(diff(sagemeans[i, ])) + 1
 if (yrs[i]==f-1) {
    offset[i] = 1
  }  else if (yrs[i]==f-2) {
    offset[i] = 2
  } else if (yrs[i] == f+1) {
    offset[i] = -1
  } else {
    offset[i] = 0
  }
}
offset[250] = -1
offset[25] = 1
ind=which(offset == 0)
ind=1:nrow(dfspat)
# plot the ts of average shrub values with proposed fire time point
fig=function(dat=sagemeans, index=dfspat$closeMatch, offs=offset, i=k, years=yrs){
  plot(dat[index[i],],type="b",main=paste0("id=",index[i],": offset=",offs[index[i]]))
  abline(v=years[index[i]]+offs[index[i]],lwd=2,col="red")
  k <<- k+1
}
k=1
fig()
####################################################################

# assemble trainign data set. use older fire in each pair as a training dataset 
n=23
ids=c(dfspat$idvar[n],dfspat$closeMatch[n])
fyear = c(dfspat$Fire_Year[ids[1]],dfspat$Fire_Year[ids[2]])


dftr = as.data.frame(sage[[ids[2]]],xy=F,na.rm = T)[,years>=fyear[2]-offset[ids[2]] ]
dftest = as.data.frame(sage[[ids[1]]],xy=F,na.rm = T)[,years>=fyear[1]-offset[ids[1]]]
matplot(t(as.matrix(dftest[tr_pred,])),type="l",main = "test")
matplot(t(as.matrix(dftr[tr,])),type="l", main="training")
set.seed(123)
tr=sample(1:nrow(dftr),round(nrow(dftr)*.02))
tr_pred=sample(1:nrow(dftest),round(nrow(dftest)*.02))

y0_train = dftr[tr,1]
y_k_prior = as.numeric(as.data.frame(sage[[ids[1]]],xy=F,na.rm = TRUE)[tr,sum(years<fyear[1]-1)])
y0_predict = dftest[tr_pred,1]
y_k_pred = as.numeric(as.data.frame(sage[[ids[2]]],xy=F,na.rm = TRUE)[tr_pred,sum(years<fyear[2]-1)])

y0_tr = ifelse(y0_train<y_k_prior,y0_train/y_k_prior,y0_train/y0_train)
y0_pred = ifelse(y0_predict<y_k_pred,y0_predict/y_k_pred,y0_predict/y0_predict)

dlist=list(y_mat=t(as.matrix(dftr[tr,])),
           N = length(tr),
           T = ncol(dftr),
           y_k_prior = y_k_prior,
           y0_tr = ifelse(is.finite(y0_tr),y0_tr,0),
           # predictions data
           y_k_pred = y_k_pred,
           y0_pred = ifelse(is.finite(y0_pred),y0_pred,0),
           N_pred = length(tr),
           T_pred = ncol(dftest),
           dftr=as.data.frame(as.matrix(dftr[tr,])),#dftr[tr,],
           dftest=dftest[tr,])


### input data
alist=dataprep(n=23,pol=dfspat,rast=sage,year=years,yroffset=offset,trsubset=.25,
               predsubset = 1,threshold = TRUE,thresholdval = .1)
matplot(t(as.matrix(alist$dftr)),type="l")
# --- fit the model to the tr fire and predict over the test fire
model2a=stan_model(paste0(path,"rdiff_discrete_landsat_non_spat_statespace.stan"))
model2=stan_model("~/Desktop/CodingClub/stan_models/rdiff_discrete_landsat_non_spat.stan")
mod2a_1=sampling(model2a, data=alist,iter=10,warmup=5,chains=1,seed=125,
             pars=c("Z","r","r_mu","r_sigma","k","k_mu","k_phi","k_sigma","eta","y_pred","y_pred01","z_pred","z_init","gamma","alpha"),
             save_warmup=FALSE)

# --- mean aboslute error
plotTraj <- function(mod, dat){
  model <- mod
  datlist <- dat
  post <- rstan::extract(model)
  yp <- apply(with(post, Z), c(2, 3), mean) #-matrix(rep(gamma,930),ncol=930)
  matplot(yp[, ], type = 'l')
}
plotTraj(mod, dat)
#matplot(t(as.matrix(alist$dftr)),type="l",add=F,col=rgb(0,0,0,.25))
plot(t(yp[-1, ]) ~ jitter(as.matrix(datlist$dftr)), pch = 3, xlab = "NLCD", ylab = "Predicted", col = rgb(0,0,0,.5))
abline(0,1,col="red",lwd=3)

y_pred <- with(post, y_pred)
err <- rep(NA, dim(y_pred)[1])

for(i in 1:length(err)){err[i] = mae(t(y_pred[i,,]),as.matrix(datlist$dftest))}
plot(density(err))
polygon(density(err), col = "blue")

### ---set up error data frame
postsample <- length(with(modlist[[1]], alpha))
err <- matrix(NA, postsample, 25)
col <- viridis(25)
for(j in 1:ncol(err)){
  post <- modlist[[j]]
  y_pred <- with(post, y_pred)
  for(i in 1:dim(y_pred)[1]){
    err[i, j] <- mae(y_pred[i, , ], as.matrix(datlist[[j]]$dftest))
  }
}
plot(density(err[, 1]), col = col[1], xlim=c(0, 9), ylim = c(0, 15), main = "Out-of-sample error", xlab="MAE")
for(i in 1:ncol(err)){polygon(density(err[, i]), col = col[i])}
#legend("topright", legend = c("simple logistic", "+varying k, +proc err", "spatial"), fill=c(col[c(1,6)],"red"))

############################################## dopar
library(foreach)
library(doParallel)
cl=makeCluster(20)
registerDoParallel(cl)
#modlist = 
foreach(i=1:25,.packages=c('raster','rstan')) %dopar% {
  #alist=dataprep(n=i,pol=dfspat,rast=sage,year=years,yroffset=offset,trsubset=.25,
  #               predsubset = 1,threshold = TRUE,thresholdval = .1)
  alist <- readRDS(paste0(path,"models_stan_eros/datasets_models_eros/dat",i,".rds"))
  temp <- stan(file = paste0(path, "rdiff_discrete_landsat_non_spat_statespace.stan"),
              data = alist, iter = 500, warmup = 450, chains = 3, seed = 125,
              pars = c("Z", "r", "r_mu", "r_sigma", "k", "k_mu", "k_phi", "k_sigma", "eta", 
                       "y_pred", "y_pred01", "z_pred", "z_init", "gamma", "alpha", "alpha_pred", "sigma_alpha_pred"),
              save_warmup=FALSE)
  saveRDS(temp, paste0(path, "models_stan_eros/mod", i, ".rds"))
}
stopCluster(cl)
#####################################################
# === explore the fit models
modlist <- list()
datlist <- list()
for(i in 1:25) {
  print(i)
  modlist[[i]] <- rstan::extract(readRDS(paste0(path, "/models_stan_eros/", "mod", i, ".rds")))
  datlist[[i]] <- readRDS(paste0(path, "/models_stan_eros/datasets_models_eros/", "dat", i, ".rds"))
}


# --- plot fire polygons with leaflet ####
poly = spTransform(df,CRS("+proj=longlat +datum=WGS84"))

r1=sage[[ids[1]]][[1]]
r2=sage[[ids[2]]][[1]]
#crs(r1) = sp::CRS("+init=epsg:3857")

pal <- colorNumeric(viridis(5), values(r1),
                    na.color = "transparent")

leaflet(poly[ids,]) %>% addPolygons(color="blue") %>%
  addTiles(group = "OSM",
                       options = providerTileOptions(minZoom = .1, maxZoom = 100)) %>%
  addProviderTiles("Esri.WorldTopoMap",    
                   group = "Topo") %>%
  addRasterImage(r1, colors = pal, opacity = 0.8,project = TRUE) %>%
  addRasterImage(r2, colors = pal, opacity = 0.8,project = TRUE) %>%
  addLegend(pal = pal, values = values(r1),
            title = "Sagebrush cover, %")



# assemble trainign data set. use older fire in each pair as a training dataset 



##### ==== set up a loop to run stan model
vec=rep(NA,25)
for(i in 1:25){
  vec[i]=modlist[[i]]$N
}
plot(vec,type="l")


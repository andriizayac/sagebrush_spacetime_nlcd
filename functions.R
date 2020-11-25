############################### This script provides functions for the workflow



# === calculate Mean Absoute Error between x and y
mae <- function(x, y){
  sum(abs(x - mean(y))) / length(x)
}

# plot the ts of average shrub values with proposed fire time point
fig <- function(dat = sagemeans, index = dfspat$idvar, offs = offset, i = k, years = yrs){
  plot(dat[i,],type="b",main=paste0("id=",index[i],": offset=",offs[index[i]], " year: ", dfspat$FireYear[i]))
  abline(v=years[i],lwd=2,col="red") #+offs[index[i]]
  k <<- k+1
}

# subset a list by logical
subsetList <- function(alist, var){
  slist = list()
  k = 1
  for(i in 1:length(alist)){
    if(var[i] == TRUE){
      slist[[k]] = alist[[i]]
    k = k+1
    }
  }
  return(slist)
}
  

# plot pre- and post- fire rasters (based on the dfspat - matched data frame) 
plotfire <- function(id = NULL, offyr = offset, focal = TRUE, pol = dfspat, rast = sage){
  par(mfrow=c(1,2))
  i = ifelse(focal==T, dfspat$idvar[id], dfspat$closeMatch[id])
  k = dfspat$Fire_Year[i] - 1984
  plot(sage[[i]][[k - 1 + offyr[i]]], main=paste0(k -1 +1984,": pre"))
  plot(sage[[i]][[k + offset[i]]], main=paste0(k + 1984,": post"))
}

# plot trajectories from a stan model
plotTraj <- function(mod, dat){
  model <- mod
  datlist <- dat
  post <- mod #rstan::extract(model)
  yp <- apply(with(post, y_pred), c(2, 3), mean) #-matrix(rep(gamma,930),ncol=930)
  matplot(yp[, ], type = 'l')
}

# prepare intput data for non-spatial stan model
dataprep <- function(n = NULL, pol = NULL, rast = NULL, years = NULL, yroffset = NULL,
                     seed = 123, trsubset = .5, predsubset = 1, threshold = FALSE, thresholdval = NULL){
  # n: #row in the matched polygon data frame
  # pol: polygon dataframe including target and ref sites
  # rast: list of rasters corresponding to each target/ref polygon (each has 33 layers=years)
  # years: vector of years - 1985:2017
  # yroffset: mismatch b/w fire year and NLCD fire; df w/ 2 cols: target, ref
  # trsubset: proportion to subset from the full dataset as training sample
  # predsubset: proportion to predict over in the test dataset
  # threshold: value [0,1] in shrub cover to include as a training dataset.Values closer to one will allow data points that introduce bias into 'r' parameter 
  
  set.seed(seed)
  n <- n
  ids <- c(pol$idvar[n], pol$closeMatch[n])
  fyear <- c(pol$Fire_Year[ids[1]], pol$Fire_Year[ids[2]])
  
  # --- the first if-else separates pairs of fires depending which one is older/younger
  if(fyear[1] < fyear[2]){
    # subset post-disturbance trajectories
    dftr <- as.data.frame(rast[[ids[1]]], xy = F, na.rm = T)[, years >= fyear[1] + yroffset[ids[1]] ]
    dftest <- as.data.frame(rast[[ids[2]]], xy = F, na.rm = T)[, years >= fyear[2] + yroffset[ids[2]] ]
    # create training indexes
    tr <- sample(1:nrow(dftr), round(nrow(dftr) * trsubset))
    tr_pred <- sample(1:nrow(dftest), round(nrow(dftest) * predsubset))
    # create training and test datasets, create pre-disutrbance means
    y0_train <- dftr[tr, 1]
    
    # --- this if-else creates a vavlue for pre-disturbance cover: previous year OR averaged since the beginning of time series
    if(sum(years < fyear[1] - 1) == 1) {
      y_k_prior <- as.numeric(as.data.frame(rast[[ids[1]]], xy = F, na.rm = T)[tr, 1:sum(years < fyear[1] - 1)])
    } else {
      y_k_prior <- as.numeric(apply(as.data.frame(rast[[ids[1]]], xy = F, na.rm = T)[tr, 1:sum(years < fyear[1] - 1)], 1, mean))
    }
    y0_predict <- dftest[tr_pred, 1]
    y_k_pred <- as.numeric(apply(as.data.frame(rast[[ids[2]]], xy = F, na.rm = T)[tr_pred, 1:sum(years < fyear[2] - 1)], 1, mean))
    
    # standardize initial values to [0,1] for non-dimensional model
    y0_tr <- ifelse(y0_train < y_k_prior, y0_train / y_k_prior, y0_train / y0_train)
    y0_pred <- ifelse(y0_predict < y_k_pred, y0_predict / y_k_pred, y0_predict / y0_predict)
    
    # --- this if-else allows to use only the pixels with initial values below a threshold 
    if(threshold == TRUE){
      y0_train <- ifelse(is.finite(y0_tr), y0_tr, 0)
      idx <- which(y0_train < thresholdval)
      # --- adds increments the threshold up in case fewer that 100 pixes fall under the threshold 
      while(length(idx) < 100) {
        thresholdval <- thresholdval + .1
        idx <- which(y0_train < thresholdval)
      }
      y_matrix <- as.matrix(dftr[tr, ])[idx, ]
      
      dlist <- list(N = nrow(y_matrix),
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

# prepare intput data for non-spatial stan model
dataprep2 <- function(n = NULL, pol = NULL, rast = NULL, years = NULL, yroffset = NULL,
                     seed = 123, trsubset = .5, predsubset = 1, threshold = FALSE, thresholdval = NULL){
  # n: #row in the matched polygon data frame
  # pol: polygon dataframe including target and ref sites
  # rast: list of rasters corresponding to each target/ref polygon (each has 33 layers = years)
  # years: vector of years - 1985:2017
  # yroffset: mismatch b/w fire year and NLCD fire; df w/ 2 cols: target, ref
  # trsubset: proportion to subset from the full dataset as training sample
  # predsubset: proportion to predict over in the test dataset
  # threshold: value [0,1] in shrub cover to include as a training dataset.Values closer to one will allow data points that introduce bias into 'r' parameter 
  
  set.seed(seed)
  n <- n
  ids <- c(pol$idvar[n], pol$closeMatch[n])
  fyear <- c(pol$Fire_Year[ids[1]], pol$Fire_Year[ids[2]])
  # --- the first if-else separates pairs of fires depending which one is older/younger
  if (pol$Fire_Year[ids[1]] > pol$Fire_Year[ids[2]]) {
    fyear <- c(pol$Fire_Year[ids[2]], pol$Fire_Year[ids[1]])
    ids <- c(pol$closeMatch[n], pol$idvar[n])
  }
 
    # subset post-disturbance trajectories
    dftr <- as.data.frame(rast[[ids[1]]], xy = F, na.rm = T)[, years >= fyear[1] + yroffset[ids[1]] ]
    dftest <- as.data.frame(rast[[ids[2]]], xy = F, na.rm = T)[, years >= fyear[2] + yroffset[ids[2]] ]
    # create training indexes
    tr <- sample(1:nrow(dftr), round(nrow(dftr) * trsubset))
    tr_pred <- sample(1:nrow(dftest), round(nrow(dftest) * predsubset))
    # create training and test datasets, create pre-disutrbance means
    y0_train <- dftr[tr, 1]
    
    # --- this if-else creates a value for pre-disturbance cover: previous year OR averaged since the beginning of time series
      if(sum(years < fyear[1] - 1) == 1) {
        y_k_prior <- as.numeric(as.data.frame(rast[[ids[1]]], xy = F, na.rm = T)[tr, 1:sum(years < fyear[1] - 1)])
      } else {
        y_k_prior <- as.numeric(apply(as.data.frame(rast[[ids[1]]], xy = F, na.rm = T)[tr, 1:sum(years < fyear[1] - 1)], 1, mean))
      }
    y0_predict <- dftest[tr_pred, 1]
    y_k_pred <- as.numeric(apply(as.data.frame(rast[[ids[2]]], xy = F, na.rm = T)[tr_pred, 1:sum(years < fyear[2] - 1)], 1, mean))
    
    # standardize initial values to [0,1] for non-dimensional model
    y0_tr <- ifelse(y0_train < y_k_prior, y0_train / y_k_prior, y0_train / y0_train)
    y0_pred <- ifelse(y0_predict < y_k_pred, y0_predict / y_k_pred, y0_predict / y0_predict)
    
    # --- this if-else allows to use only the pixels with initial values below a threshold (whithin equilibria \in [0, 1])
    if(threshold == TRUE){
      y0_train <- ifelse(is.finite(y0_tr), y0_tr, 0)
      idx <- which(y0_train < thresholdval)
      # --- adds increments the threshold up in case fewer that 100 pixes fall under the threshold 
        while(length(idx) < 100) {
          thresholdval <- thresholdval + .1
          idx <- which(y0_train < thresholdval)
        }
      y_matrix <- as.matrix(dftr[tr, ])[idx, ]
      
      dlist <- list(N = nrow(y_matrix),
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
  return(dlist)
}

# prepare input data for non-spatial stan model without predictions
dataprep3 <- function(alist, k, df, tr_prop = .5, thresholdval = 1){
  rast = alist[[k]]
  fyear = dfspat$FireYear[k]
  # subset post-disturbance trajectories
  dftr <- as.data.frame(rast, xy = F, na.rm = T)[, years >= fyear]
  # create training indexes
  tr <- sample(1:nrow(dftr), round(nrow(dftr) * tr_prop))
  # create training and test datasets, create pre-disutrbance means
  y0_train <- dftr[tr, 1]
  y0_pred <- dftr[-tr, 1]
  
  # --- this if-else creates a value for pre-disturbance cover: previous year OR averaged since the beginning of time series
  if(sum(years < fyear - 1) == 1) {
    y_k_prior <- as.numeric(as.data.frame(rast, xy = F, na.rm = T)[tr, 1:sum(years < fyear)])
    y_k_pred <- as.numeric(as.data.frame(rast, xy = F, na.rm = T)[-tr, 1:sum(years < fyear)])
  } else {
    y_k_prior <- as.numeric(apply(as.data.frame(rast, xy = F, na.rm = T)[tr, 1:sum(years < fyear[1])], 1, mean))
    y_k_pred <- as.numeric(apply(as.data.frame(rast, xy = F, na.rm = T)[-tr, 1:sum(years < fyear[1])], 1, mean))
  }
  
  # standardize initial values to [0,1] for non-dimensional model. scales either to a prop of x/k or returns 1 
  y0_tr <- ifelse(y0_train < y_k_prior, y0_train / y_k_prior, y0_train / y0_train)
  y0_pr <- ifelse(y0_pred < y_k_pred, y0_pred / y_k_pred, y0_pred / y0_pred)

  # --- this part allows to use only the pixels with initial values below a threshold (whithin equilibria \in [0, 1])
    y0_train <- ifelse(is.finite(y0_tr), y0_tr, 0)
    y0_pred <- ifelse(is.finite(y0_pr), y0_pr, 0)
    idx <- which(y0_train < thresholdval)
    # --- adds increments the threshold up in case fewer that 100 pixes fall under the threshold 
    while(length(idx) < 100) {
      thresholdval <- thresholdval + .1
      idx <- which(y0_train < thresholdval)
    }
    y_matrix <- as.matrix(dftr[tr, ])[idx, ]
    
    dlist <- list(N = nrow(y_matrix),
                  tr = tr,
                  T = ncol(dftr),
                  # training data
                  y_mat=t(y_matrix),
                  y_k_prior = y_k_prior[idx],
                  y0 = y0_train[idx],
                  # prediction data
                  n = nrow(dftr[-tr, ]),
                  y0pred = y0_pred,
                  y_k_pred = y_k_pred,
                  # record threshold used
                  cover_threshold = thresholdval,
                  # input data sets
                  dftr=as.data.frame(y_matrix))

  return(dlist)
}


i = 1
alist1 <- dataprep3(sage,1,dfspat,.15,.25)
alist1 <- dataprep(n = i, pol = dfspat, rast = sage, year = years, yroffset = offset, seed = 123,
                  trsubset = .1,
                  predsubset = 1,threshold = TRUE,thresholdval = .1)
alist2 <- dataprep2(n = i, pol = dfspat, rast = sage, year = years, yroffset = offset, seed = 123, 
                  trsubset = .1,
                  predsubset = 1,threshold = TRUE,thresholdval = .1)

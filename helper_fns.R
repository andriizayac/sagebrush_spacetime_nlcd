# === change color alpha
# credit: R. McElreath
col_a <- function(col, alpha) {
  c <- col2rgb(col)
  cola <- rgb(c[1]/255, c[2]/255, c[3]/255, alpha)
  return(cola)
}

# === mask functions for sf
st_erase = function(x, y) {
  st_difference(
    st_geometry(x) %>% st_buffer(0), 
    st_union(st_combine(st_geometry(y))) %>% st_buffer(0)
  )
}
st_eraseq = function(x, y) { 
  st_difference(x, st_union(st_combine(st_geometry(y))) )
}

# === vectorize function from: 'matrixcalc'
vec <- function(x) {return( t(t(as.numeric(x))))}

# === Data extracting functions
glm.dat <- function(tsage, tfires, tpxlcov, i=NULL, clN = NULL) {
  # this function generates data input for log-log gompertz Poisson, Gamma or Gaussian model
  # input: sage cover rasters, fire dataset, pixel-level covaraites, polygons index, number of clusters to be used
  # output: response - sage % for t = 2,3,..T
  #   predictor - negative sage %, encoding 0 -> .25, for t = 1,2,...T-1
  #   group - cl - pixel-level classes for random effect
  mat <- data.matrix(tsage[[i]][,c(tfires$FireYer[i]-1984):31])
  matcl <- outer(tpxlcov[[i]][, paste0("cluster", clN)], rep(1, ncol(mat)))
  xpr <- vec(mat[,-ncol(mat)])
  ypr <- vec(mat[,-1])
  outs <- which(xpr == 0)
  if(length(outs) == 0) {
    cl <- vec(matcl[,-1])
    dat <- data.frame(y = ypr,
                      x = xpr,
                      cl = as.factor(cl))
  } else {
    cl <- vec(matcl[,-1])[-outs]
    dat <- data.frame(y = ypr[-outs],
                      x = xpr[-outs],
                      cl = as.factor(cl))
  }
  return(dat)
}
glm.dat.init <- function(tsage, tfires, tpxlcov, i=NULL, clN = NULL){
  mat <- tsage[[i]][,c(tfires$FireYer[i]-1984):31]
  T = ncol(mat)
  colnames(mat) <- 1:T-1
  mat$cl = as.factor(tpxlcov[[i]][,paste0("cluster", clN)])
  dat.init <- pivot_longer(mat, cols = c(1:T)) %>%
    rename(t = name, cover = value) %>% 
    mutate(t = as.numeric(t)) %>%
    filter(t <= 5)
  return(dat.init)
}

# === Gomperz analytical solution in y(t):
# y(t) = K * exp(C * exp(-alpha*t)) 
gomp <- function(a, b, n0, t) {
  # vector-valued Gompertz model
  # in: coefs from log-log and N0-Poisson models and a vector of time steps
  # out: matrix of values nrow = length(a), and ncol = length(t)
  
  K <- exp(-a/b)
  C <- log(n0/K)
  u <- K * exp(C * exp( outer(b, t) ))
  return(u)
}

# === calculate pixel Mahalanobis distance
pxl.dist <- function(dat, dat.m, coefs, k) {
  # finds a match for each focal pixel to a centroid of the matching cluster in the reference site
  # dat: focal pixel-level covariates
  # dat.m: reference pixel-level covariates
  # k: number of cluster to be considered
  # output: an integer vector for each pixel in a target dataset
  clusterv <- paste0("cluster", k) 
  df <- dat.m %>% 
    dplyr::select(1:5, clusterv) %>%  
    group_by_at(clusterv) %>%
    summarize_all(mean) %>% 
    # mutate(ninit = coefs[[j]]$n0) %>% 
    dplyr::select(-starts_with("cluster"), -stab) 
  df0 <- dat %>% 
    dplyr::select(1:5, clusterv) %>% 
    rename(cl = names(.)[6]) %>% 
    # left_join(coefs[[i]][, 3:4], by = "cl") %>% 
    # mutate_at(vars(n0),~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)) %>% 
    # rename(ninit = n0) %>% 
    # dplyr::select(-ninit) %>% 
    group_by(cl) %>% 
    summarize_all(mean) %>% dplyr::select(-cl, -stab)
  
  scdf <- scale(rbind(df0, df))
  n <- nrow(scdf)
  mah.p <- mah.dist(scdf[1:(n-k), ], scdf[(n-k+1):n, ])

  pdist <- apply(mah.p, 1, which.min)
  
  return(pdist)
}

# === Error functions
mae <- function(e, rel = FALSE, k = NULL){
  # inputs (vector-valued)
    # yhat: predicted values
    # y: observed values
    # rel: error relative or not
    # k: a value to relativize by
  # outputs
    # a scalar error value
  
  if(rel == TRUE) {
    sum( abs(e) / k) / length(e)  
  } else {
    sum(abs(e)) / length(e)
  }
}
rmse <- function(e, rel = FALSE, k = k){
  # inputs (vector-valued)
  # yhat: predicted values
  # y: observed values
  # rel: error relative or not
  # k: a value to relativize by
  # outputs
  # a scalar error value
  
  if(rel == TRUE) {
      sqrt( sum( ((e)^2) / k) / length(e))   
  } else {
    sqrt( sum( (e)^2 ) / length(e) )
  }
}
bias <- function(e) {
  sum(e) / length(e)
}
# --- end error functions

# === modified mahalanobis distance from StatMatch with decreased tolerance
# --- full credit to StatMatch::mahalanobis.dist
# --- not used
mah.dist <- function(data.x, data.y = NULL, vc = NULL) 
{
  xx <- as.matrix(data.x)
  if (is.null(data.y)) 
    yy <- as.matrix(data.x)
  else yy <- as.matrix(data.y)
  if (is.null(vc)) {
    if (is.null(data.y)) 
      vc <- var(xx)
    else vc <- var(rbind(xx, yy))
  }
  ny <- nrow(yy)
  md <- matrix(0, nrow(xx), ny)
  for (i in 1:ny) {
    md[, i] <- mahalanobis(xx, yy[i, ], cov = vc, tol = 1e-19)
  }
  if (is.null(data.y)) 
    dimnames(md) <- list(rownames(data.x), rownames(data.x))
  else dimnames(md) <- list(rownames(data.x), rownames(data.y))
  sqrt(md)
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

# === prep data for the NULL (ie regional) model
# --- combine data sets - different number of years post-fire don't matter because model is fit recursively
dat.gen.null.full <- function(datin, tfires, i) {
  # this function differs from the one used in the glm in that zeros from the predictor are not thrown away and are stored in the data
  mat <- nullsage[[i]][,c(tfires$FireYer[i]-1984):33]
  T <- ncol(mat)
  colnames(mat) <- 1:T-1
  dat <- mat %>% 
    pivot_longer(cols = 1:T) %>%
    rename(t = name, cover = value) %>% 
    mutate(t = as.numeric(t)) 
  df <- data.frame(y = dat$cover[dat$t != 0],
                   x = dat$cover[dat$t != max(dat$t)],
                   t = dat$t[dat$t != max(dat$t)], 
                   id = i)
  return(df)
}

dat.gen.null <- function(datin, tfires, i) {
  mat <- nullsage[[i]][,c(tfires$FireYer[i]-1984):33]
  T <- ncol(mat)
  colnames(mat) <- 1:T-1
  dat <- mat %>% 
    pivot_longer(cols = 1:T) %>%
    rename(t = name, cover = value) %>% 
    mutate(t = as.numeric(t)) 
  df <- data.frame(y = dat$cover[dat$t != 0],
                   x = dat$cover[dat$t != max(dat$t)],
                   t = dat$t[dat$t != max(dat$t)], 
                   id = i) %>% 
    filter(x != 0)
  return(df)
}

dat.gen.null.init <- function(datin, tfires, i=NULL){
  mat <-datin[[i]][,c(tfires$FireYer[i]-1983):33]
  T = ncol(mat)
  colnames(mat) <- 1:T-1
  dat <- pivot_longer(mat, cols = c(1:T)) %>%
    rename(t = name, cover = value) %>% 
    mutate(t = as.numeric(t)) %>% 
    filter(t <= 5) %>% 
    mutate(id = i)
  return(dat)
}

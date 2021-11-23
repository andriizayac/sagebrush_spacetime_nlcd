# === mask function for sf
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

# === Data extracting function for linear model
lm.dat <- function(tsage, tfires, tpxlcov, i=NULL, clN = NULL){
  mat <- tsage[[i]][,c(tfires$FireYer[i]-1984):31]
  T = ncol(mat)
  colnames(mat) <- 1:ncol(mat)-1
  mat$cl = as.factor(tpxlcov[[i]][,paste0("cluster", clN)])
  dat.init <- pivot_longer(mat, cols = c(1:T)) %>%
    rename(t = name, cover = value) %>% 
    mutate(t = as.numeric(t)) 
  return(dat.init)
}

# === Gomperz equation
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

# === calculate pixel Mahalanobis distance
pxl.dist <- function(dat, dat.m, coefs, k) {
  # finds a match for each focal pixel to a centroid of the matching cluster in the reference site
  # dat: focal pixel-level covariates
  # dat.m: reference pixel-level covariates
  # k: number of cluster to be considered
  # output: an integer vector for each pixel in a target dataset
  clusterv <- paste0("cluster", k) 
  df <- dat.m %>% dplyr::select(1:3, clusterv) %>%  
    group_by_at(clusterv) %>%
    summarize_all(mean) %>% 
    mutate(ninit = coefs[[j]]$n0) %>% 
    dplyr::select(-starts_with("cluster"), -ninit) 
  df0 <- dat %>% 
    dplyr::select(1:3, clusterv) %>% 
    rename(cl = names(.)[4]) %>% 
    left_join(coefs[[i]][, 3:4], by = "cl") %>% 
    mutate_at(vars(n0),~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)) %>% 
    rename(ninit = n0) %>% 
    dplyr::select(-cl, -ninit) 
  
  scdf <- scale(rbind(df0, df))
  n <- nrow(scdf)
  mah.p <- StatMatch::mahalanobis.dist(scdf[1:(n-k), ], scdf[(n-k+1):n, ])
  
  pdist <- apply(mah.p, 1, which.min)
  
  return(pdist)
}
pxl.dist.m <- function(dat, dat.m, k) {
  # finds a match for each focal pixel to a centroid of the matching cluster in the reference site
  # dat: focal pixel-level covariates
  # dat.m: reference pixel-level covariates
  # k: number of cluster to be considered
  # output: an integer vector for each pixel in a target dataset
  clusterv <- paste0("cluster", k) 
  df <- dat.m %>% map(. %>% 
                        select(1:5, clusterv) %>% 
                        group_by_at(clusterv) %>% 
                        summarize_all(mean) %>% 
                        dplyr::select(-starts_with("cluster")))
  
  df1 <- do.call(rbind, lapply(seq_along(df), function(x) {
    df[[x]]$X1 <- x
    df[[x]]
  })) 
  mah.p <- StatMatch::mahalanobis.dist(dat[, 1:5], df1[,1:5])
  
  pdist <- apply(mah.p, 1, which.min)
  
  return(pdist)
}
pxl.dist.wt <- function(dat, dat.m, k) {
  # finds a set of matching clusters in a ref and calculates weights for pixels within min + sd of mah distance
  # dat: focal pixel-level covariates
  # dat.m: reference pixel-level covariates
  clusterv <- paste0("cluster", k) 
  df <- dat.m %>% dplyr::select(-starts_with("cluster"))
  
  ind <- 1:nrow(df)
  if(nrow(df) > 10000) {
    set.seed(123)
    ind <- sample(ind, 10000)
    df <- df[ind,]
  }
  mah.p <- StatMatch::mahalanobis.dist(dat[, 1:5], df)
  
  pdist <- apply(mah.p, 1, function(x) {
    dat.m[ind,clusterv][which(x <= min(x) + sd(x) )]
  } )
  
  weights <- sapply(pdist, function(x) table(x)/length(x))
  clusters <- sapply(weights, function(x) as.numeric(names(x)) )
  
  return( list(clusters = clusters, weights = weights) )
}



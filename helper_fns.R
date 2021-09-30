
# === Data extracting functions
glm.dat <- function(tsage, tfires, tpxlcov, i=NULL, clN = NULL) {
  # this function generates a data input for log-log gompertz Poisson or Gaussian model
  # input: sage cover rasters, fire dataset, pixel-level covaraites, polygons index, number of clusters to be used
  # output: response - sage % for t = 2,3,..T
  #   predictor - negative sage %, encoding 0 -> .25, for t = 1,2,...T-1
  #   group - cl - pixel-level classes for random effect
  mat <- data.matrix(tsage[[i]][,c(tfires$FireYer[i]-1984):33])
  xpr <- ifelse(vec(mat[,-ncol(mat)]) == 0, .2, vec(mat[,-ncol(mat)]))
  ypr <- ifelse(vec(mat[,-1]) == 0, .2, vec(mat[,-1]))
  dat <- data.frame(y = log(ypr/xpr),
                    x = log(xpr),
                    cl = as.factor(rep(tpxlcov[[i]][, paste0("cluster", clN)], times = ncol(mat)-1)))
  return(dat)
}
glm.dat.init <- function(tsage, tfires, tpxlcov, i=NULL, clN = NULL){
  mat <- tsage[[i]][,c(tfires$FireYer[i]-1984):33]
  T = ncol(mat)
  colnames(mat) <- 1:T-1
  mat$cl = as.factor(tpxlcov[[i]][,paste0("cluster", clN)])
  dat <- pivot_longer(mat, cols = c(1:T)) %>%
    rename(t = name, cover = value) %>% 
    mutate(t = as.numeric(t), 
           cover = ifelse(cover == 0, runif(1e6, 0,1), cover)) %>%
    filter(t <= 5)
  return(dat)
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
pxl.dist <- function(dat, dat.m, k) {
  # finds a match for each focal pixel to a cluster in the reference site
  # - replaces cluster value for matched cluster
  # dat: focal pixel-level covariates
  # dat.m: reference pixel-level covariates
  clusterv <- paste0("cluster", k) 
  df <- dat.m %>% dplyr::select(-starts_with("cluster"))
    # group_by_at(clusterv) %>%
    # summarize_all(mean) %>% dplyr::select(-starts_with("cluster"))
    ind <- 1:nrow(df)
    if(nrow(df) > 10000) {
      ind <- sample(ind, 10000)
      df <- df[ind,]
    }
  mah.p <- StatMatch::mahalanobis.dist(dat[, 1:5], df)
  
  pdist <- apply(mah.p, 1, function(x) {
    dat.m[,clusterv][which(x <= min(x) + sd(x) )]
    } )
  
  weights <- sapply(pdist, function(x) table(x)/length(x))
  clusters <- sapply(weights, function(x) as.numeric(names(x)) )
  
  # dat$cluster <- dat.m[ind, clusterv][ apply(mah.p, 1, which.min) ]
  return( list(clusters = clusters, weights = weights) )
}



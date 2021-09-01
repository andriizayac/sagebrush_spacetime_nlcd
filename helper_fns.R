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

# This script calculates MAE(%), RMSE(%), Bias based on the output from `model_pred.R`
pkgs <- c("tidyverse")
sapply(pkgs, require, character.only = T)

source("helper_fns.R")
# === load data
tpxlcov <- readRDS("data/tpxlcov.rds")
tdfEnv <- readRDS("data/tdfEnv_covars.rds") %>% 
        dplyr::select(tmax, tmin, dem_mean, chili_mean, sm03mean, sm03sd, ppt)
tfires <- readRDS("data/tfires.rds")
tsage <- readRDS("data/tsage.rds")

N <- nrow(tfires)
k <- sapply(tpxlcov, function(x) { mean(x$prefire) } )
# --- end

# === set up the error matrix
# rows: 1-MAE; 2-RMSE; 3-MAE%; 4-RMSE%; 5-SD
dferr <- data.frame(matrix(NA, nr = 5, nc = 15))
names(dferr) <- c("Null", "Site", paste(2:14))

# === Null predictions
errnull <- read.csv("outputs/null_dat_pred_T9.csv")
errs <- readRDS("outputs/matchpred.rds")
y <- sapply(errs$obs, function(x) { x[[10]] } )
en <- errnull$yhat - y

dferr$Null[1] <- mae(en, FALSE)
dferr$Null[2] <- rmse(en, FALSE)
dferr$Null[3] <- mae(en, rel = TRUE, k)
dferr$Null[4] <- rmse(en, rel = TRUE, k)
dferr$Null[5] <- sd(en)
# --- end null
        
# === Site level predictions
errs <- readRDS("outputs/matchpred.rds")
y <- sapply(errs$obs, function(x) { x[[10]] } )
yhat <- sapply(errs$pred, function(x) { x[[10]] } )
es <- yhat - y

dferr$Site[1] <- mae(es, FALSE)
dferr$Site[2] <- rmse(es, FALSE)
dferr$Site[3] <- mae(es, rel = TRUE, k)
dferr$Site[4] <- rmse(es, rel = TRUE, k)
dferr$Site[5] <- sd(es)

# --- end site

# === Site + cluster level predictions
errsc <- readRDS("outputs/matchsitecluster.rds")
ys = lapply(errsc, function(x) {
        lapply(x$obs, function(y) { y[,10] } )} )
yhats = lapply(errsc, function(x) {
                lapply(x$pred, function(y) { y[,10] } )} )
escmat <- matrix(NA, nc = 13, nr = N)
# escmat <- array(NA, c(13, N, 14) )
for(M in 3:14){
        ly = ys[[M]]
        lyhat = yhats[[M]]
        for(i in 1:N){
                escmat[i, M-1] = mean(lyhat[[i]] - ly[[i]], na.rm = TRUE)
        }
}

for(M in 3:14){ 
        # ins = which(is.finite(escmat[,M-1])) # for in-sample predictions with bad models
        ins = 1:N 
        dferr[1, M+1] = mae(escmat[ins,M-1], FALSE)
        dferr[2, M+1] = rmse(escmat[ins,M-1], FALSE)
        dferr[3, M+1] = mae(escmat[ins,M-1], rel = TRUE, k[ins])
        dferr[4, M+1] = rmse(escmat[ins,M-1], rel = TRUE, k[ins])
        dferr[5, M+1] = sd(escmat[ins,M-1])
}
matplot(t(as.matrix(dferr[1:4,])), type = "l")

efin <- data.frame(Null = en, Site = es) %>% 
        bind_cols(as.data.frame(escmat[,-1]))
write.csv(efin, "outputs/error_bar_N_df.csv", row.names = FALSE)

dferr %>%
        mutate(across(where(is.numeric), round, 3)) %>%
        mutate(Metric = c("MAE", "RMSE","MAE%", "RMSE%", "SD")) %>%
        relocate(Metric, Null) %>%
        write.csv(., file = "figures/Table2.csv", row.names = FALSE)

# put tables together
dat.out <- read.csv("figures/Table2b.csv")
dat.in <- read.csv("figures/Table2a.csv") 
dat.prop <- dat.out
for(i in 2:ncol(dat.prop)){dat.prop[,i] = dat.in[,i] - dat.out[,i]}
# write.csv(dat.prop, file = "figures/Table2c.csv", row.names = FALSE)



# === BIAS ===
hn <- read.csv("outputs/null_dat_pred.csv")
hs <- readRDS("outputs/matchpred.rds")
hsc <- readRDS("outputs/matchsitecluster.rds")

# --- observed and predicted nlcd time series
ymat <- matrix(NA, nr = N, nc = 29)
yerrmat <- array(NA, c(15, N, 29))
ins <- lengths(hs$obs)
for(i in 1:ncol(ymat)){
        idx = which(ins >= i)
        ymat[idx,i] = sapply(hs$obs[idx], function(x) { x[[i]] } )
}
for(i in 1:ncol(ymat)){
        idx = which(ins >= i)
        # --- null and site levels
        yerrmat[1,idx,i] = hn$yhat[hn$t == i-1][idx] - ymat[idx,i]
        yerrmat[2,idx,i] = sapply(hs$pred[idx], function(x) {x[[i]]} ) - ymat[idx,i]
        # --- cluster levels
        ys = lapply(hsc, function(x) {
                        lapply(x$obs[idx], function(y) { y[,i] } )} )
        yhats = lapply(hsc, function(x) {
                        lapply(x$pred[idx], function(y) { y[,i] } )} )
        for(M in 3:length(hsc)){
                ly = ys[[M]]
                lyhat = yhats[[M]]
                for(n in 1:length(idx)) {
                        yerrmat[M+1,idx[n],i] = mean(lyhat[[n]] - ly[[n]], na.rm = TRUE)
                }
        }
}

saveRDS(list(ymat = ymat, yerrmat = yerrmat), file = "outputs/bias.rds")


# === Recovery to pre-wildfire state 10 years after a wildfire
errs <- readRDS("outputs/matchpred.rds")
y <- sapply(errs$obs, function(x) { x[[10]] } )
yhat <- sapply(errs$pred, function(x) { x[[10]] } )
prek <- k
recov0 <- y/k
recov <- yhat/k

cat("NLCD observed recovery at 10 years: ", mean(recov0), sd(recov0))
cat("Predicted recovery at 10 years: ", mean(recov), sd(recov))

cat("Predicted proportion P(Y_hat > k) 10 years after disturbance: ", 
    sum(yhat > k)/length(k))


pckg <- c("ggplot2", "dplyr", "tidyr")
sapply(pckg, require, character.only = T)

N <- length(nullout$mae)

# === Null model: requires loading outputs/null_model_lm.rds and outputs/null_model.rds files
par(mfrow = c(1, 2))
plot(density(mae), lwd = 2, xlab = "Mean absolute error, %", main = "")
plot(mae ~ jitter(tsf), pch = 19, col = rgb(0, 0, 0, .25), 
     ylab = "Mean absolute error, %", xlab = "Time since fire, years")
mtext("Null Model", side = 3, line = -2, outer = TRUE)
dev.off()
# ---
data.frame(Gompertz = mae, LinearFit = maelm) %>% 
  pivot_longer(cols = 1:2) %>% 
  rename(MAE = value, model = name) %>% 
  ggplot(aes(MAE, group = model, fill = model)) + geom_density(adjust=1.5, alpha=.5) +
  labs(y = "Density") +
  theme_bw()
# --- global models vs pre-disturbance 
plot(mae ~ kvec, pch = 19, col = rgb(.5,0,0,.75), bty = "n", xlab = "Pre-disturbance mean, %", ylab = "MAE")
points(maelm ~ kvec, pch = 19, col = rgb(0,.5,0,.75) )
abline(0, 1, lwd = 2)

# --- null MAPE and MAE 
nullout <- readRDS("outputs/null_model.rds")
tsf <- nullout$tsf
mapenull <- rep(NA, N)

for(i in 1:N) {
  y <- nullout$datnull[[i]]
  yhat <- nullout$yhatnull[[i]]
  mapenull[i] = mean( abs(y[, tsf[i]] - yhat[, tsf[i]])/y[, tsf[i]])
}
# --- null MAE
maenull <- rep(NA, N)
for(i in 1:N) {
  y <- nullout$datnull[[i]]
  yhat <- nullout$yhatnull[[i]]
  maenull[i] = mean( abs(y[, tsf[i]] - yhat[, tsf[i]]) )
}

par(mfrow = c(1, 2), mar = c(4, 4, 4, 1))
plot(density(mapenull), bty = "n", lwd = 3, xlab = "MAPE", main = "")
plot(density(maenull), bty = "n", lwd = 3, xlab = "MAE", main = "")
mtext("Null Model", side = 3, line = -2, outer = TRUE)


# === 





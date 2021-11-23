# === plot cases
lout <- list()
for(i in 1:N){
  lout[[i]] <- data.matrix(tsage[[i]][,c(tfires$FireYer[i]-1984):33])
}
k <- 1
f <- function(k) {
  matplot(t(lout[[k]]), type = "l", main = k)
  k <<- k + 1
}
f(k)

# === fit brms
i = 1
dat1 <- glm.dat(tsage, tfires, tpxlcov, i, 4)
glm1 <- glmer(y ~ (1|cl) + (0+log(x)|cl), 
              offset =log(x), family = poisson(), data = dat1) #

j = 200
dat2 <- glm.dat(tsage, tfires, tpxlcov, j, 4)
glm2 <- glmer(y ~ (1|cl) + (0+x|cl), family = poisson(), data = dat2) #

fit2 <- dat2 %>% 
  #mutate(x = log(x)) %>% 
  update(fit1, newdata = .)

# === plot raw data
mat1 <- data.matrix(tsage[[i]][,c(tfires$FireYer[i]-1984):31])
mat2 <- data.matrix(tsage[[j]][,c(tfires$FireYer[j]-1984):31])
matplot(t(mat1), type = "l", col = rgb(0,0,0,.5))
matplot(t(mat2), type = "l")
# --- add empirical means
cl <- tpxlcov[[i]]$cluster4
df <- mat1 %>% 
  as.data.frame() %>% 
  mutate(cl = cl) %>% 
  group_by(cl) %>% 
  summarize_all(mean) %>% 
  dplyr::select(-cl)
apply(df, 1, lines, lwd = 2, col = "magenta")
# add predictions
gpred <- gomp(coefs[[i]]$a, coefs[[i]]$b, coefs[[i]]$n0, 0:c(ncol(df) - 1))
for(h in 1:nrow(gpred)) { lines(gpred[h,], lwd = 3, col = "purple")}
# --- 
fe <- summary(fit1)$fixed[,1]
u = gomp(fe[1], fe[2], exp(-1), 0:30)[1,]
lines(u, lwd = 3)

re <- ranef(fit1)$cl[,1,]

ure <- gomp(fe[1] + re[,1], fe[2] + re[,2], exp(-1), 0:12)
ure <- gomp(fe + re[,1], re[,2], exp(-1), 0:30)
apply(ure, 1, lines, lwd = 2, col = "brown")



dat1 %>% 
  mutate(x = -1*log(x), y = log(y)) %>% 
  ggplot(aes(x = x, y = y)) + geom_point()

# === brms models
dat1 <- glm.dat(tsage, tfires, tpxlcov, 40, 4)
fit1 <- dat1 %>% 
  mutate(nlogx = -log(x)) %>% # 
  brm(y ~ 1 + log(x) + offset(log(x)),
      data = ., family = poisson(), 
      #prior = c(prior(normal(0, 1), class = "Intercept"),
      #          prior(normal(0, .5), class = "b", lb = -1, ub = 0)),
      iter = 2000, chains = 3, cores = 3 #,sample_prior = "only"
  )

fit2 <- dat1 %>% 
  mutate(nlogx = -log(x)) %>% # 
  brm(y ~ 0 + Intercept + log(x) + offset(log(x)),
      data = ., family = poisson(), 
      #prior = c(prior(normal(0, 1), class = "Intercept"),
      #          prior(normal(0, .5), class = "b", lb = -1, ub = 0)),
      iter = 2000, chains = 3, cores = 3 
  )

fit3 <- dat1 %>% 
  mutate(nlogx = -log(x)) %>% # 
  brm(y ~ 1 + log(x) + (1 + log(x)|cl) + offset(log(x)),
      data = ., family = poisson(), 
      #prior = c(prior(normal(0, 1), class = "Intercept"),
      #          prior(normal(0, .5), class = "b", lb = -1, ub = 0)),
      iter = 2000, chains = 3, cores = 3 
  )

fit4 <- dat1 %>% 
  mutate(nlogx = -log(x)) %>% # 
  brm(y ~ 0 + Intercept + log(x) + (Intercept + log(x)|cl) + offset(log(x)),
      data = ., family = poisson(), 
      #prior = c(prior(normal(0, 1), class = "Intercept"),
      #          prior(normal(0, .5), class = "b", lb = -1, ub = 0)),
      iter = 2000, chains = 3, cores = 3 
  )

fit5 <- dat1 %>% 
  mutate(nlogx = -log(x)) %>% # 
  brm(y ~ 0 + (Intercept + log(x)|cl) + offset(log(x)),
      data = ., family = poisson(), 
      #prior = c(prior(normal(0, 1), class = "Intercept"),
      #          prior(normal(0, .5), class = "b", lb = -1, ub = 0)),
      iter = 2000, chains = 3, cores = 3 
  )

fit6 <- dat1 %>% 
  mutate(nlogx = -log(x)) %>% # 
  brm(y ~ 0 + (1 + log(x)|cl) + offset(log(x)),
      data = ., family = poisson(), 
      #prior = c(prior(normal(0, 1), class = "Intercept"),
      #          prior(normal(0, .5), class = "b", lb = -1, ub = 0)),
      iter = 2000, chains = 3, cores = 3 
  )


# =======
fit1npneg <- dat1 %>% 
  mutate(nlogx = -log(x)) %>% # 
  brm(y ~ 1 + nlogx + offset(-log(x)),
      data = ., family = poisson(), 
      #prior = c(prior(normal(0, 1), class = "Intercept"),
      #          prior(normal(0, .5), class = "b", lb = -1, ub = 0)),
      iter = 100, chains = 2, cores = 2 #,sample_prior = "only"
  )

fit1pr0 <- dat1 %>% 
  mutate(nlogx = -log(x)) %>% # 
  brm(y ~ 1 + nlogx + offset(-log(x)),
      data = ., family = poisson(), 
      prior = c(prior(exponential(10), class = "Intercept"),
                prior(lognormal(2, 3), class = "b", lb = 0.01)),
      iter = 1000, chains = 2, cores = 2 #,sample_prior = "only"
  )

fit1pr1 <- dat1 %>% 
  mutate(nlogx = -log(x)) %>% update(fit1pr1, newdata = ., iter = 1000, chains = 2, cores = 2)
  brm(y ~ 1 + log(x) + offset(-log(x)),
      data = ., family = poisson(), 
      prior = c(prior(exponential(10), class = "Intercept")),
                #prior(lognormal(2, 3), class = "b")),
      iter = 100, chains = 2, cores = 2 #,sample_prior = "only"
  )

# ------------
  i = 283
dat <- glm.dat(tsage, tfires, tpxlcov, i, 4)
ffix <- brmsformula(y ~ 1 + x)
fre <- brmsformula(y ~ 1 + x + (1 + x|cl))
  
f <- brm(ffix, family = poisson, data = dat, chains = 2, cores = 2, iter = 200)
  
datu <- glm.dat(tsage, tfires, tpxlcov, 40, 4)
fu <- update(f, newdata = dat, fre)


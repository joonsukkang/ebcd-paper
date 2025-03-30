library(here)
library(tidyverse)
source(here('code', 'utils', 'fit_utils.R'))

draw.fig <- TRUE


# V0 and Sigma
V0 <- matrix(0, nrow=500, ncol=2)
V0[  1:10, 1] <- 1/sqrt(10)
V0[ 11:20, 2] <- 1/sqrt(10)

p <- nrow(V0)
K <- ncol(V0)

Sigma <- V0 %*% diag(x=c(399, 299)) %*% t(V0) + diag(x=rep(1, p))


# generate and save samples
n <- 50
nseeds <- 50

mat <- lapply(1:nseeds, function(x){
  set.seed(x)
  X <- MASS::mvrnorm(n = n, mu = rep(0, ncol(Sigma)), Sigma = Sigma)
  return(X)
})

mat <- do.call(rbind, mat)
MASS::write.matrix(mat, 
                   file = here('data', 'sim', 'sim1.csv'), 
                   sep = ',')
rm(mat)


# fit models
mat <- read_csv(file = here('data', 'sim', 'sim1.csv'), 
                col_names=FALSE, show_col_types = FALSE)
mat <- as.matrix(mat); colnames(mat) <- NULL


t0 <- Sys.time()
res.pca    <- run.on.mat(mat, function(x){ fit.pca(x,K)$L }, n, p, nseeds)
t1 <- Sys.time()
print((t1 - t0)/nseeds)

t0 <- Sys.time()
res.ebcd_pl <- run.on.mat(mat, function(x){ fit.ebcd_pl(x,K)$L }, n, p, nseeds)
t1 <- Sys.time()
print((t1 - t0)/nseeds)

t0 <- Sys.time()
res.ebcd_l <- run.on.mat(mat, function(x){ fit.ebcd_l(x,K)$L }, n, p, nseeds)
t1 <- Sys.time()
print((t1 - t0)/nseeds)

t0 <- Sys.time()
res.ebmf_pl <- run.on.mat(mat, function(x){ fit.ebmf_pl(x,K)$L }, n, p, nseeds)
t1 <- Sys.time()
print((t1 - t0)/nseeds)

t0 <- Sys.time()
res.ebmf_l <- run.on.mat(mat, function(x){ fit.ebmf_l(x,K)$L }, n, p, nseeds)
t1 <- Sys.time()
print((t1 - t0)/nseeds)

t0 <- Sys.time()
res.l1ppca <- run.on.mat(mat, function(x){ fit.l1ppca(x,K)$L }, n, p, nseeds)
t1 <- Sys.time()
print((t1 - t0)/nseeds)

t0 <- Sys.time()
res.spc    <- run.on.mat(mat, function(x){ fit.spc(x,K)$L }, n, p, nseeds)
t1 <- Sys.time()
print((t1 - t0)/nseeds)



### ebpca
# run run_ebpca_sim1.py
res.ebpca <- as.matrix(read_csv(here('output', 'sim', 'ebpca_sim1.csv'), 
                                col_names = FALSE, show_col_types = FALSE))
### gpower, choose parameter by smallest dcov
# run run_gpower_sim1.m
res.gpower <- opt.gpower(simno = 1, n, nseeds, K, V0, Sigma, draw.fig)

 
# compute distance
df.sim1 <- data.frame()
methods <- c('pca', 'gpower', 'ebcd_pl', 'ebcd_l', 'ebmf_pl', 'ebmf_l', 'l1ppca', 'spc', 'ebpca')
for (method in methods){
  df.sim1 <- rbind(df.sim1, 
               data.frame(method = method, 
                          seed = 1:nseeds, 
                          dist.type = rep(c(paste0('d', 1:K), 'dOR', 'dcov'), each = nseeds),
                          dist = c(dist.on.res(get(paste0('res.', method))/sqrt(n), V0, Sigma, nseeds, K))))
}

saveRDS(df.sim1, file=here('output', 'sim', 'df.sim1.rds'))

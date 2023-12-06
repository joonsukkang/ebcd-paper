library(here)
library(tidyverse)
source(here('code', 'fit_utils.R'))

draw.fig <- TRUE


# V0 and Sigma
V0 <- matrix(0, nrow=500, ncol=3)
V0[  1:10, 1] <- 1/sqrt(10)
V0[ 11:50, 2] <- 1/sqrt(40)
V0[ 51:150, 3] <- 1/sqrt(100)


p <- nrow(V0)
K <- ncol(V0)

Sigma <- V0 %*% diag(x=c(9, 7, 4)) %*% t(V0) + diag(x=rep(1, p))


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
                   file = here('data', 'sim2.csv'), 
                   sep = ',')
rm(mat)


# fit models
mat <- read_csv(file = here('data', 'sim2.csv'), 
                col_names=FALSE, show_col_types = FALSE)
mat <- as.matrix(mat); colnames(mat) <- NULL


t0 <- Sys.time()
res.pca    <- run.on.mat(mat, function(x){ fit.pca(x,K)$L }, n, p, nseeds)
t1 <- Sys.time()
print((t1 - t0)/nseeds)

t0 <- Sys.time()
res.ebcd <- run.on.mat(mat, function(x){ fit.ebcd(x,K)$L }, n, p, nseeds)
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
system("cd '/Users/jkang/Library/CloudStorage/Box-Box/research/ebcd-paper/code/';
       /Users/jkang/miniconda3/bin/python3 run_ebpca_sim2.py",
       ignore.stdout=FALSE, wait=TRUE)
res.ebpca <- as.matrix(read_csv(here('output', 'ebpca_sim2.csv'), 
                                col_names = FALSE, show_col_types = FALSE))
# gpower: hard to specify mu

# compute distance
df.sim2 <- data.frame()
methods <- c('pca', 'ebcd', 'l1ppca', 'spc', 'ebpca')
for (method in methods){
  df.sim2 <- rbind(df.sim2, 
               data.frame(method = method, 
                          seed = 1:nseeds, 
                          dist.type = rep(c(paste0('d', 1:K), 'dOR', 'dcov'), each = nseeds),
                          dist = c(dist.on.res(get(paste0('res.', method))/sqrt(n), V0, Sigma, nseeds, K))))
}

saveRDS(df.sim2, file=here('output', 'df.sim2.rds'))

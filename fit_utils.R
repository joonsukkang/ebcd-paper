# functions to fit data
source(here('code', 'L1PPCA.R'))
library(PMA)
library(ebcd)


fit.pca <- function(X, K){
  svdX <- svd(X, nu=K, nv=K)
  Z <- svd(X)$u[,1:K, drop=FALSE]
  L <- svd(X)$v[,1:K, drop=FALSE] %*% diag(x=svd(X)$d[1:K], nrow=K)
  out.list <- list(Z=Z, L=L)
  
  return(out.list)
}

fit.irlba <- function(X, K){
  fit <- irlba::irlba(X, nv=K, nu=0)
  L <- fit$v %*% diag(x=fit$d[1:K], nrow=K)
  out.list <- list(L=L)
  
  return(out.list)
}


fit.ebcd <- function(X, K){

  fit <- ebcd(X=X, Kmax=K)
  out.list <- list(Z=fit$Z/sqrt(nrow(X)), L=fit$EL*sqrt(nrow(X)))

  return(out.list)
}


fit.l1ppca <- function(X, 
                       K,
                       lambda.grid = seq(0, 10, len=25)){
  

  lambda.eq <- matrix(rep(lambda.grid, times=K), ncol=K)
  fit.cv <- l1ppca.cv(X, lambda.mat = lambda.eq)
  if(max(fit.cv$opt.lambda) == max(lambda.grid)){
    warning('maximum lambda chosen as optimal; expand the grid')
  }
  fit <- l1ppca(X, lambda = fit.cv$opt.lambda)
  
  out.list <- list(Z=fit$Z, L=fit$L)
  
  return(out.list) 
}


ortho <- function(U, seed=0){
  
  set.seed(seed)
  
  mat <- matrix(rnorm(nrow(U)^2), nrow=nrow(U))
  mat[, 1:ncol(U)] <- U
  mat <- qr.Q(qr(mat))
  ortho <- mat[,-(1:ncol(U))]
  return(ortho)
}


fit.spc <- function(X, K, sumabsvs = seq(1, sqrt(ncol(X)), len = 25)){
  
  
  d <- c()
  u <- matrix(nrow = nrow(X), ncol = 0)
  v <- matrix(nrow = ncol(X), ncol = 0)
  
  for (k in 1:K){
    if(k==1){
      X.target <- X
      cv.out <- SPC.cv(x = X.target, 
                       sumabsvs = sumabsvs, trace = FALSE, center = FALSE)
      out <- SPC(x = X.target, 
                 sumabsv = cv.out$bestsumabsv, K = 1, trace = FALSE, center = FALSE)
      d <- out$d
      u <- matrix(out$u, ncol=1)
      v <- matrix(out$v, ncol=1)
    } else{
      X.target <- t(ortho(u)) %*% (X - u %*% diag(x=d, ncol=length(d)) %*% t(v))
      cv.out <- SPC.cv(x = X.target, 
                       sumabsvs = sumabsvs, trace = FALSE, center = FALSE)
      out <- SPC(x = X.target, 
                 sumabsv = cv.out$bestsumabsv, K = 1, trace = FALSE, center = FALSE)
      d <- c(d, out$d)
      u <- cbind(u, ortho(u) %*% matrix(out$u, ncol=1))
      v <- cbind(v, out$v)
    }
  }
  
  out.list <- list(Z = u, L = v %*% diag(d, nrow=K))
  
  return(out.list)
}


# function to run method over all data sets
run.on.mat <- function(mat, method1, n, p, nseeds){
  matr <- matrix(0, nrow = 0, ncol = p)
  
  for (seed in 1:nseeds){
    
    row_from <- (seed-1) * n + 1
    row_to <- seed * n
    mat1 <- mat[row_from:row_to, ]
    
    matr <- rbind(matr, t(method1(mat1))) 
  }
  return(matr)
}


# function to compute distance
calc.dist <- function(L.by.sqrtN, V0 = NULL, Sigma){
  
  d.vec <- c()
  if(!is.null(V0)){
    
    # distance by vector
    L.remaining <- apply(L.by.sqrtN, 2, function(x){
      if(sum(x^2)>0){x/sqrt(sum(x^2))}else x
    } )
    for (k in 1:ncol(L.by.sqrtN)){
      d <- acos(abs(c(t(V0[,k] %*% L.remaining))))/(pi/2)
      d.vec <- c(d.vec, min(d))
      if( k < ncol(V0) ){
        L.remaining <- L.remaining[, -which.min(d)]
      }
    }
    
    # subspace difference
    L <- qr.Q(qr(L.by.sqrtN)) # convert to orthogonal basis before computing distance
    V0 <- qr.Q(qr(V0))
    d.vec <- c(d.vec, 
               norm(t(polar(t(L) %*% V0)) %*% t(V0) - t(L), 'f'))
  }
  
  # covariance difference
  d.vec <- c(d.vec, 
             norm(Sigma - tcrossprod(L.by.sqrtN), 'f'))
  
  return(d.vec)
}

# compute distance over all results
dist.on.res <- function(res, V0 = NULL, Sigma, nseeds, K){
  
  
  if(!is.null(V0)){
    matr <- matrix(0, nrow = 0, ncol = K+2)
  }else{
    matr <- matrix(0, nrow = 0, ncol = 1)
  }
  
  
  for (seed in 1:nseeds){
    
    row_from <- (seed-1) * K + 1
    row_to <- seed * K
    res1 <- res[row_from:row_to, ,drop=FALSE]
    
    matr <- rbind(matr, calc.dist(L.by.sqrtN = t(res1), V0, Sigma))
  }
  if(!is.null(V0)){
    colnames(matr) <- c(paste0('d',1:K), 'dOR', 'dcov')
  }
  
  return(matr)
}


opt.gpower <- function(simno, n, nseeds, K, V0, Sigma, draw.fig=FALSE){
  res.gpower <- read_csv(here('output', paste0('gpower_sim', simno, '.csv')), 
                         col_names = FALSE, show_col_types = FALSE)
  res.gpower <- as.matrix(res.gpower)/sqrt(n)
  dist.gpower <- data.frame()
  for (g in 1:99){
    g_from <- ((g-1)*(K*nseeds)+1)
    g_to <- g*(K*nseeds)
    matr.g <- res.gpower[g_from:g_to, ]
    dist.g <- colMeans(dist.on.res(matr.g, V0, Sigma, nseeds, K))
    dist.gpower <- rbind(dist.gpower,
                         data.frame(g=g, t(dist.g)))
  }
  if(draw.fig){
    dist.gpower %>%
      ggplot()+
      geom_line(aes(x=g, y=d1, col='d1'))+
      geom_line(aes(x=g, y=d2, col='d2'))+
      geom_line(aes(x=g, y=dOR, col='dOR'))+
      geom_line(aes(x=g, y=dcov, col='dcov'))+
      geom_vline(aes(xintercept=which.min(dcov)), col='black')+
      ggtitle(paste0('GPower: simulation ', simno))
  }
  print(dist.gpower[which.min(dist.gpower$dcov), ])
  
  g <- which.min(dist.gpower$dcov)
  g_from <- ((g-1)*(K*nseeds)+1)
  g_to <- g*(K*nseeds)
  res.gpower.opt  <- res.gpower[g_from:g_to, ]

  return(res.gpower.opt)
}


calc.mse <- function(X.test, L){
  
  L <- L[, !is.nan(colSums(L^2))]
  Lortho <- qr.Q(qr(L))
  X.proj <- X.test %*% tcrossprod(Lortho)
  mse <- mean((X.test - X.proj)^2)
  return(mse)
}

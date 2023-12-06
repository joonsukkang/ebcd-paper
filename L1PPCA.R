library(irlba)
polar <- function(A){
  svdA <- svd(A)
  out <- svdA$u %*% t(svdA$v)
  return(out)
}

soft.thresholding <- function(a, lambda){
  out <- sign(a) * pmax(0, abs(a) - lambda)
  return(out)
}

l1ppca <- function(X, 
                   lambda,
                   maxiter = 500, 
                   tol=1e-5,
                   init.Z=NULL){

  K <- length(lambda)  
  if(is.null(init.Z)){
    svdX <- irlba(X, nu=K, nv=0)
    Z <- svdX$u
  }else{
    Z <- init.Z
  }
  
  L <- matrix(0, nrow=ncol(X), ncol=K)
  obj.old <- Inf
  vec.obj <- c()
  
  for (iter in 1:maxiter){
    
    # Shrinkage step
    L <- t(X) %*% Z
    for (k in 1:K){
      L[,k] <- soft.thresholding(L[,k], lambda[[k]])
    }
    
    # Rotation step
    Z <- polar(X %*% L) 
    
    # check convergence
    obj <- sum((X - Z %*% t(L))^2) + 2*sum(lambda * colSums(abs(L)))
    vec.obj <- c(vec.obj, obj)
    if(obj - obj.old > - tol) break
    obj.old <- obj
  }
  
  out.list <- list(Z=Z, L=L, vec.obj=vec.obj)
  return(out.list)
}


l1ppca.cv <- function(X,
                      lambda.mat, # each row is a candidate value vector
                      maxiter = 500,
                      tol = 1e-5,
                      nfold = 5,
                      seed = 0){
  
  set.seed(seed)
  row.idx <- sample(nrow(X))
  mat.mse <- matrix(0, nrow = nrow(lambda.mat), ncol = nfold)
  
  for (fold in 1:nfold){
    idx.from <- (fold - 1) * floor(nrow(X)/nfold) + 1
    idx.to   <- fold * floor(nrow(X)/nfold)
    
    X.train <- X[row.idx[-(idx.from:idx.to)], ]
    X.test <- X[row.idx[idx.from:idx.to], ]
    
    svdX <- irlba(X.train, nu=ncol(lambda.mat), nv=0)
    init.Z <- svdX$u
    
    for (i in 1:nrow(lambda.mat)){
      
      fit <- l1ppca(X.train, lambda.mat[i,], maxiter = maxiter, tol = tol,
                    init.Z = init.Z)
      Lortho <- qr.Q(qr(fit$L))
      X.proj <- (X.test %*% Lortho) %*% t(Lortho)
      mse <- mean((X.test - X.proj)^2)
      mat.mse[i, fold] <- mse
    }
  }
  
  vec.mse <- rowMeans(mat.mse)
  opt.lambda <- lambda.mat[which.min(vec.mse),]
  out.list <- list(opt.lambda = opt.lambda, 
                   vec.mse = vec.mse, 
                   mat.mse = mat.mse)
  return(out.list)
}




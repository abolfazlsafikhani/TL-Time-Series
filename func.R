#-------------------------------------------------------------------------------
library(MASS)
library(CVXR)
library(glmnet)
library(plotrix)

###### FISTA low-rank plus sparse estimating functions ######
fista.LpS <- function(A, b, lambda, mu, group, K, niter, backtracking = TRUE, x.true = NULL){
  #A      : design matrix. matrix consisting of all N samples in rows
  #b      : predictor matrix. 
  #lambda : tunning parameter: penalty for sparse matrix
  #mu     : tunning parameter: penalty for lowrank matrix
  #group  : a vector of starting point for each segment
  #K      : number of group
  #niter  : maximum of iteration
  tnew = t <- 1
  p <- dim(A)[2]; N <- dim(A)[1]
  xnew2 = y2 <- matrix(0, nrow = p, ncol = p)
  ls_y1 <- list(); ls_x1 <- list()
  for(i in 1:K){
    ls_y1[[i]] <- matrix(0, nrow = p, ncol = p)
    ls_x1[[i]] <- matrix(0, nrow = p, ncol = p)
  }
  ls_xnew1 <- ls_x1
  ls_AtA <- list()
  ls_Atb <- list()
  group <- c(group, N+1)
  for(i in 1:K){
    ls_AtA[[i]] <- t(A[group[i]:group[i+1]-1,]) %*% A[group[i]:group[i+1]-1,]  
    ls_Atb[[i]] <- t(A[group[i]:group[i+1]-1,]) %*% b[group[i]:group[i+1]-1,]
  }
  Atb <- t(A) %*% b
  
  if(backtracking == TRUE){
    L <- norm(A, "F")^2 / 5
    gamma <- 2
  }else{
    L <- norm(A, "F")^2
  }
  
  obj.val = rel.err <- c()
  ls_prox1 <- list()
  for(i in 1:niter){
    if(backtracking == TRUE){
      L.bar <- L
      found <- FALSE
      while(found == FALSE){
        #for(l in 1:200){
        for(i in 1:K){
          y <- ls_y1[[i]] + y2
          ls_prox1[[i]] <- prox.sparse.func1(ls_y1[[i]], y, A[group[i]:group[i+1]-1,], b[group[i]:group[i+1]-1,], 2*L.bar, lambda, ls_AtA[[i]], ls_Atb[[i]]) 
        }
        
        prox2 <- prox.nuclear.func1(y2, ls_y1, A, b, 2*L.bar, mu, ls_AtA, ls_Atb, K)
        
        #print(prox2)
        ### Restricted solution space
        for(j in 1:p){
          for(k in 1:p){
            if(abs(prox2[j,k]) > 0.25){
              prox2[j,k] <- 0.25 * sign(prox2[j,k])
            }
          }
        }
        
        #prox <- prox1 + prox2
        if(f.func1(ls_prox1, prox2, A, b, group, K) <= Q.func1(ls_prox1, prox2, ls_y1, y2, A, b, L.bar, ls_AtA, ls_Atb, group, K)){
          found <- TRUE
        }else{
          L.bar <- L.bar * gamma
        }
        #L.bar <- L.bar * gamma
      }
      L <- L.bar
    }
    ls_x1 <- ls_xnew1 
    x2 <- xnew2
    ls_xnew1 <- ls_prox1
    xnew2 <- prox2
    t = tnew
    tnew <- (1 + sqrt(1 + 4*t^2))/2
    for(i in 1: K){
      ls_y1[[i]] <- ls_xnew1[[i]] + (t - 1) / tnew * (ls_xnew1[[i]] - ls_x1[[i]])  
    }
    y2 <- xnew2 + (t - 1) / tnew * (xnew2 - x2)
    #xnew <- xnew1 + xnew2
    
    #obj.val <- c(obj.val, f.func(xnew, A, b) + sparse.pen(xnew1, lambda) + nuclear.pen(xnew2, mu))
    #rel.err <- c(rel.err, norm(xnew - x.true, "F") / norm(x.true, "F"))
  }
  #return(list(sparse.comp = xnew1, lr.comp = xnew2, obj.val = obj.val, rel.err = rel.err))
  return(list(sparse.comp = ls_xnew1, lr.comp = xnew2))
}

fista.Lasso<- function(A, b, lambda, mu, group, K, niter, backtracking = TRUE, x.true = NULL){
  tnew = t <- 1
  p <- dim(A)[2]; N <- dim(A)[1]
  xnew2 = y2 <- matrix(0, nrow = p, ncol = p)
  ls_y1 <- list(); ls_x1 <- list()
  for(i in 1:K){
    ls_y1[[i]] <- matrix(0, nrow = p, ncol = p)
    ls_x1[[i]] <- matrix(0, nrow = p, ncol = p)
  }
  ls_xnew1 <- ls_x1
  ls_AtA <- list()
  ls_Atb <- list()
  group <- c(group, N+1)
  for(i in 1:K){
    ls_AtA[[i]] <- t(A[group[i]:group[i+1]-1,]) %*% A[group[i]:group[i+1]-1,]  
    ls_Atb[[i]] <- t(A[group[i]:group[i+1]-1,]) %*% b[group[i]:group[i+1]-1,]
  }
  Atb <- t(A) %*% b
  
  if(backtracking == TRUE){
    L <- norm(A, "F")^2 / 5
    gamma <- 2
  }else{
    L <- norm(A, "F")^2
  }
  
  obj.val = rel.err <- c()
  ls_prox1 <- list()
  for(i in 1:niter){
    if(backtracking == TRUE){
      L.bar <- L
      found <- FALSE
      while(found == FALSE){
        #for(l in 1:200){
        for(i in 1:K){
          y <- ls_y1[[i]] #+ y2
          ls_prox1[[i]] <- prox.sparse.func1(ls_y1[[i]], y, A[group[i]:group[i+1]-1,], b[group[i]:group[i+1]-1,], 2*L.bar, lambda, ls_AtA[[i]], ls_Atb[[i]]) 
        }
        
        #prox <- prox1 + prox2
        if(f.func1(ls_prox1, matrix(0,p,p), A, b, group, K) <= Q.func1(ls_prox1, matrix(0,p,p), ls_y1, matrix(0,p,p), A, b, L.bar, ls_AtA, ls_Atb, group, K)){
          found <- TRUE
        }else{
          L.bar <- L.bar * gamma
        }
        #L.bar <- L.bar * gamma
      }
      L <- L.bar
    }
    ls_x1 <- ls_xnew1 
    #x2 <- xnew2
    ls_xnew1 <- ls_prox1
    #xnew2 <- prox2
    t = tnew
    tnew <- (1 + sqrt(1 + 4*t^2))/2
    for(i in 1: K){
      ls_y1[[i]] <- ls_xnew1[[i]] + (t - 1) / tnew * (ls_xnew1[[i]] - ls_x1[[i]])  
    }
    #y2 <- xnew2 + (t - 1) / tnew * (xnew2 - x2)
    #xnew <- xnew1 + xnew2
    
    #obj.val <- c(obj.val, f.func(xnew, A, b) + sparse.pen(xnew1, lambda) + nuclear.pen(xnew2, mu))
    #rel.err <- c(rel.err, norm(xnew - x.true, "F") / norm(x.true, "F"))
  }
  #return(list(sparse.comp = xnew1, lr.comp = xnew2, obj.val = obj.val, rel.err = rel.err))
  return(list(sparse.comp = ls_xnew1, lr.comp = xnew2))
}


fista.nuclear <- function(A, b, lambda, mu, group, K, niter, backtracking = TRUE, x.true = NULL){
  tnew = t <- 1
  p <- dim(A)[2]; N <- dim(A)[1]
  xnew2 = y2 <- matrix(0, nrow = p, ncol = p)
  ls_y1 <- list(); ls_x1 <- list()
  ls_prox1 <- list()
  for(i in 1:K){
    ls_y1[[i]] <- matrix(0, nrow = p, ncol = p)
    ls_x1[[i]] <- matrix(0, nrow = p, ncol = p)
    ls_prox1[[i]] <- matrix(0, nrow = p, ncol = p)
  }
  #ls_xnew1 <- ls_x1
  ls_AtA <- list()
  ls_Atb <- list()
  group <- c(group, N+1)
  for(i in 1:K){
    ls_AtA[[i]] <- t(A[group[i]:group[i+1]-1,]) %*% A[group[i]:group[i+1]-1,]  
    ls_Atb[[i]] <- t(A[group[i]:group[i+1]-1,]) %*% b[group[i]:group[i+1]-1,]
  }
  Atb <- t(A) %*% b
  
  if(backtracking == TRUE){
    L <- norm(A, "F")^2 / 5
    gamma <- 2
  }else{
    L <- norm(A, "F")^2
  }
  
  obj.val = rel.err <- c()
  
  for(i in 1:niter){
    if(backtracking == TRUE){
      L.bar <- L
      found <- FALSE
      while(found == FALSE){
        #for(l in 1:200){
        #for(i in 1:K){
        #  y <- ls_y1[[i]] + y2
        #  ls_prox1[[i]] <- prox.sparse.func1(ls_y1[[i]], y, A[group[i]:group[i+1]-1,], b[group[i]:group[i+1]-1,], 2*L.bar, lambda, ls_AtA[[i]], ls_Atb[[i]]) 
        #}
        
        prox2 <- prox.nuclear.func1(y2, ls_y1, A, b, 2*L.bar, mu, ls_AtA, ls_Atb, K)
        
        #print(prox2)
        ### Restricted solution space
        for(j in 1:p){
          for(k in 1:p){
            if(abs(prox2[j,k]) > 0.25){
              prox2[j,k] <- 0.25 * sign(prox2[j,k])
            }
          }
        }
        
        #prox <- prox1 + prox2
        if(f.func1(ls_prox1, prox2, A, b, group, K) <= Q.func1(ls_prox1, prox2, ls_y1, y2, A, b, L.bar, ls_AtA, ls_Atb, group, K)){
          found <- TRUE
        }else{
          L.bar <- L.bar * gamma
        }
        #L.bar <- L.bar * gamma
      }
      L <- L.bar
    }
    #ls_x1 <- ls_xnew1 
    x2 <- xnew2
    #ls_xnew1 <- ls_prox1
    xnew2 <- prox2
    t = tnew
    tnew <- (1 + sqrt(1 + 4*t^2))/2
    #for(i in 1: K){
    #  ls_y1[[i]] <- ls_xnew1[[i]] + (t - 1) / tnew * (ls_xnew1[[i]] - ls_x1[[i]])  
    #}
    y2 <- xnew2 + (t - 1) / tnew * (xnew2 - x2)
    #xnew <- xnew1 + xnew2
    
    #obj.val <- c(obj.val, f.func(xnew, A, b) + sparse.pen(xnew1, lambda) + nuclear.pen(xnew2, mu))
    #rel.err <- c(rel.err, norm(xnew - x.true, "F") / norm(x.true, "F"))
  }
  #return(list(sparse.comp = xnew1, lr.comp = xnew2, obj.val = obj.val, rel.err = rel.err))
  return(list(lr.comp = xnew2))
}


shrinkage <- function(y, tau){
  z <- matrix(0, nrow = nrow(y), ncol = ncol(y))
  for(i in 1:nrow(y)){
    for(j in 1:ncol(y)){
      z[i,j] <- sign(y[i,j]) * max(abs(y[i,j]) - tau, 0)
    }
  }
  return(z)
}
shrinkage.lr <- function(y, tau){
  z <- rep(0, length(y))
  for(i in 1:length(y)){
    z[i] <- sign(y[i]) * max(0, abs(y[i]) - tau)
  }
  return(z)
}
gradf.func <- function(x, AtA, Atb){
  return(AtA %*% x - Atb)
}
nuclear.pen <- function(x, lambda){
  d <- svd(x)$d
  return(lambda * sum(d))
}
sparse.pen <- function(x, lambda){
  return(lambda*sum(x))
}
f.func1 <- function(ls_prox1, prox2, A, b, group, K){
  SUM <- 0
  for(i in 1:K){
    SUM <- SUM + (0.5 * norm(A[group[i]:group[i+1]-1,] %*% (ls_prox1[[i]] + prox2) - b[group[i]:group[i+1]-1,], "F")^2)  
  }
  return(SUM)
}
Q.func1 <- function(ls_prox1, prox2, ls_y1, y2, A, b, L, ls_AtA, ls_Atb, group, K){
  SUM1 <- f.func1(ls_y1, y2, A, b, group, K)
  SUM2 <- 0
  for (i in 1:K){
    SUM2 <- SUM2 + sum((ls_prox1[[i]] + prox2 - ls_y1[[i]] - y2)*gradf.func(ls_y1[[i]] + y2, ls_AtA[[i]], ls_Atb[[i]]))
  }
  SUM3 <- 0
  for(i in 1:K){
    SUM3 <- 0.5 * L * norm(ls_prox1[[i]] + prox2 - ls_y1[[i]] - y2, "F")^2  
  }
  return(SUM1 + SUM2 + SUM3)
}
prox.nuclear.func1 <- function(w1, ls_y, A, b, L, lambda, ls_AtA, ls_Atb, K){
  C <- w1
  
  for (i in 1:K){
    C <- C - (1 / L) * gradf.func(ls_y[[i]] + w1, ls_AtA[[i]], ls_Atb[[i]])
  }
  Y <- C
  d <- shrinkage.lr(svd(Y)$d, 2*lambda / L)
  return(svd(Y)$u %*% diag(d) %*% t(svd(Y)$v))
}
prox.sparse.func1 <- function(w1, y, A, b, L, lambda, AtA, Atb){
  Y <- w1 - (1 / L) * gradf.func(y, AtA, Atb)
  
  return(shrinkage(Y, 2*lambda / L))
}
obj.func <- function(x.lr, x.sparse, A, b, lambda, mu){
  ### x.sparse is a list
  m <- length(x.sparse)
  loss <- 0
  for(i in 1:m){
    loss <- loss + f.func((x.lr[[i]] + x.sparse[[i]]), A, b) + sparse.pen(x.sparse[[i]], lambda) + nuclear.pen(x.lr[[i]], mu)
  }
  return(loss)
}


#lasso_debais: For calcualating entry-wise confidence interval.
lasso_debais <- function(X, Y, Sigma_e, alpha=2, beta, p, d, mu = NULL){
  # Compute confidence intervals and p-values.
  #
  # Args:
  #   X        :  matrix consisting of all N samples in rows (dimension: N*pd)
  #   Y        :  matrix consisting of all N samples in rows (dimension: N*p)
  #   Sigma_e  :  covariance matrix for Sigma(dimension: p*p)
  #   alpha    :  episodes growth rate
  #   beta     :  estimated coefficient matrix(dimension: p*pd)
  #   p        :  dimension for time-series
  #   d        :  VAR(d) model for time-series
  #   mu       :  tunning paramter. NULL for finding by itself(takes a lot of time). mu=1~2 is an empircal choice 
  #
  # Returns:
  #   V: conditional variance in our paper
  #   beta: center point of confidence interval
  #confidence intercal: [beta - qnorm(0.975)*V/sqrt(n), beta + qnorm(0.975)*V/sqrt(n)]
  
  n <- nrow(X)
  n1 <- floor(sqrt(n))
  n_grid <- c(n1)
  n_cur = n1
  i = 1
  alpha = 2
  while(floor(n_cur+alpha^i) < n){
    n_grid = c(n_grid,floor(alpha^i))
    n_cur = floor(n_cur +alpha^i)
    i = i+1
  }
  n_grid = c(n_grid,n-n_cur)
  k = i-1
  
  Sigma <- matrix(0,p*d,p*d)
  start = 0
  E = diag(1,p*d)
  Mat_ls = list(matrix(0,p*d,p*d))
  for(i in 1:(k+1)){
    n_cur = n_grid[i]
    n_all = sum(n_grid[1:i])
    Sigma = (t(X[(start+1):(start+n_cur),])%*%X[(start+1):(start+n_cur),]+Sigma*(n_all-n_cur))/n_all
    start = start+n_cur
    Mat_l <- c()
    for(j in 1:(p*d)){
      if(is.null(mu)){
        m=Variable(p*d);ea=rep(0,p*d);ea[j]=1
        mu_l = 2*(solve(Problem(Minimize(p_norm(ea-Sigma%*%m,Inf))), solver="SCS")$value)
      }else{
        mu_l = sqrt(2*log(p*d)/n_all)/mu
      }
      Mat_l <-cbind(Mat_l,optimfunc(Sigma, j, mu_l, maxiter=50, threshold=1e-2,  m_start=Mat_ls[[i]][j,])$optsol)
    }
    Mat_ls = c(Mat_ls,list(t(Mat_l)))
  }
  
  E_es <- (Y - X%*%beta)
  beta_debais <- beta
  for(i in 1:(k+1)){
    seg_cur <- (sum(n_grid[1:i])+1):sum(n_grid[1:(i+1)])
    if(length(seg_cur) > 1){
      beta_debais = beta_debais + Mat_ls[[i+1]]%*%t(X[seg_cur,])%*%E_es[seg_cur,]/n  
    }else{
      beta_debais = beta_debais + Mat_ls[[i+1]]%*%t(matrix(X[seg_cur,],nrow = 1))%*%matrix(E_es[seg_cur,], nrow = 1)/n 
    }
    
  }
  
  V <- rep(0,p*d)
  for(i in 1:(k+1)){
    seg_cur <- (sum(n_grid[1:i])+1):sum(n_grid[1:(i+1)])
    if(length(seg_cur) > 1){
      mat_cur = Mat_ls[[i+1]]%*%t(X[seg_cur,])  
    }else{
      mat_cur = Mat_ls[[i+1]]%*%t(matrix(X[seg_cur,], nrow = 1))
    }
    V = V + apply(mat_cur,1,square_sum)/n
  }
  V_all <- c()
  for(i in 1:p){
    V_all = cbind(V_all,V*Sigma_e[i,i])
  }
  output = list(V = V_all, beta = beta_debais)
  output
}


#calaulate FPR,TPR,coverage rate
#false positive rate: V/U
#true positive rate: S/T
POS_rate <- function(fit, beta_es , beta, n){
  V = 0
  S = 0
  U = 0
  T = 0
  cover = 0
  len = qnorm(0.975)*sqrt(fit$V)/sqrt(n)
  LB = fit$beta - len
  UB = fit$beta + len
  for(i in 1:nrow(beta)){
    for(j in 1:ncol(beta)){
      if(beta[i,j]<UB[i,j] &&beta[i,j]>LB[i,j]){
        cover = cover + 1
      }
      if(beta[i,j]==0){
        U = U+1
        if(0>UB[i,j] || 0<LB[i,j]){
          V = V+1
        }
      }else{
        T = T+1
        if(0>UB[i,j] || 0<LB[i,j]){
          S = S+1
          #print(c(i,j))
        }
      }
    }
  }
  output = c(V/U,S/T,cover/(U+T))
}


#VAR Oracle Trans-lasso
#n1: target sample size; n2: auxiliary sample size
tran_learn <- function(X, Y, lambda1,lambda2, p, d, n1, n2){
  beta <- c()
  #firststep:
  for(i in 1:p){
    #lambda1 = sqrt(2*log(p*d)/n2)
    Y_temp = Y[,i]
    la_temp = glmnet(X, Y_temp, lambda=lambda2,family = "gaussian", intercept = F, alpha=1) 
    beta <- cbind(beta, la_temp$beta)
  }
  #second step
  Y0 = Y[1:n1,]
  X0 = X[1:n1,]
  Z0 = Y0 - X0%*%beta
  delta <- c()
  for(i in 1:p){
    #lambda2 = sqrt(2*log(p*d)/n1)
    Z_temp = Z0[,i]
    la_temp = glmnet(X0, Z_temp, lambda=lambda1,family = "gaussian", intercept = F, alpha=1) 
    delta <- cbind(delta, la_temp$beta)
  }
  beta_tran <- beta + delta
  beta_tran
}



lasso <- function(X , Y , p , d, lambda){
  beta <- c()
  for(i in 1:p){
    #lambda1 = sqrt(2*log(p*d)/(N1-tran))
    Y_temp = Y[,i]
    la_temp = glmnet(X, Y_temp, lambda=lambda,family = "gaussian", intercept = F, alpha=1) 
    beta <- cbind(beta, la_temp$beta)
  }
  beta
}

lasso.cv <- function(X, Y, p, d, lambda.grid, fold = 5){
  l <- length(lambda.grid)
  cv_error <- rep(0,l)
  
  for(j in 1:fold){
    X_new <- c(); Y_new <- c(); index_test <- c(0, 0)
    X_test <- c(); Y_test <- c()
    
    start <- 1; end <- nrow(X)
    segment <- floor((end - start)/fold)
    index_test[1] <- (j-1)*segment+1; index_test[2] <- j*segment+1
    index_fit <- (start:end)[-(index_test[1] : index_test[2])]
    index_test <- (start:end)[(index_test[1] : index_test[2])] 
    
    X_new <- rbind(X_new, X[index_fit,]); Y_new <- rbind(Y_new, Y[index_fit,])
    X_test <- rbind(X_test, X[index_test,]); Y_test <- rbind(Y_test, Y[index_test,])
    
    for(i in 1:l){
      beta <- lasso(X_new , Y_new , p , d, lambda.grid[i])
      cv_error[i] <- cv_error[i] + sum((Y_test - X_test%*%beta)^2)
    }
  }
  
  index <- which.min(cv_error)
  beta <- lasso(X , Y , p , d, lambda.grid[index])
  beta
}


#m0 -- start point, L -- radius
optimfunc <- function ( sigma, i, mu, maxiter=50, threshold=1e-2, 
                        m_start=NULL) {
  p <- nrow(sigma);
  if(is.null(m_start)){m_start=rep(0,p) }
  beta=m_start
  diff.norm2 <- 1;
  last.norm2 <- 1;
  iter <- 1;
  iter.old <- 1;
  beta.old <- beta;
  sigma.tilde <- sigma;
  diag(sigma.tilde) <- 0;
  vs <- -sigma.tilde%*%beta; 
  
  while ((iter <= maxiter) && (diff.norm2 >= threshold*last.norm2)){    
    
    for (j in 1:p){
      oldval <- beta[j];
      v <- vs[j];
      if (j==i)
        v <- v+1;    
      beta[j] <- SoftThreshold(v,mu)/sigma[j,j];
      if (oldval != beta[j]){
        vs <- vs + (oldval-beta[j])*sigma.tilde[,j];
      }
    }
    
    iter <- iter + 1;
    if (iter==2*iter.old){
      d <- beta - beta.old;
      diff.norm2 <- sqrt(sum(d*d));
      last.norm2 <-sqrt(sum(beta*beta));
      iter.old <- iter;
      beta.old <- beta;
      if (iter>10)
        vs <- -sigma.tilde%*%beta;
    }
  }
  
  returnlist <- list("optsol" = beta, "iter" = iter)
  return(returnlist)
}


SoftThreshold <- function( x, lambda ) {
  # Standard soft thresholding
  if (x>lambda){
    return (x-lambda);}
  else {
    if (x< (-lambda)){
      return (x+lambda);}
    else {
      return (0); }
  }
}

VAR_tau<-function(A,Sigma_e,p){
  matrix(solve(diag(1,p^2) - kronecker(A,A))%*%matrix(c(Sigma_e),ncol=1),p,p)
}

square_sum <- function(x){
  sum(x**2)
}

#a special case of generate_sigma function(sigma = diag(1,p)).
generate <- function(N,p,d,A){
  epi <- matrix(rnorm(N*p),N,p) 
  X <- matrix(0,N,p)
  for(i in 1:d){
    X[i,] <- epi[i,]
  }
  for(i in (d+1):N){
    for(j in 1:d){
      X[i,] <- A[,((j-1)*p+1):(j*p)]%*%X[i-j,] + X[i,]
    }
    X[i,] <- X[i,]+epi[i,]
  }
  X
}


#generate time series according to Sgima(covariance matrix of error term)
generate_sigma <- function(N,p,d,A,sigma, initX = NULL){
  epi <- mvrnorm(n = N+d, rep(0, p), sigma)
  #epi <- matrix(runif(N*p,-2,2),N,p) 
  X <- matrix(0,N+d,p)
  if(is.null(initX)){
    for(i in 1:d){
      X[i,] <- epi[i,]
    }
  }else{
    for(i in 1:d){
      X[i,] <- initX[(i-1)*p+1,i*p]
    }
  }
  
  for(i in (d+1):(N+d)){
    for(j in 1:d){
      X[i,] <- A[,((j-1)*p+1):(j*p)]%*%X[i-j,] + X[i,]
    }
    X[i,] <- X[i,]+epi[i,]
  }
  X[(d+1):(N+d),]
}

#generate time series according to epi(error term)
generate_epi <- function(N,p,d,A,epi, initX = NULL){#epi -- N+d * p 
  #epi <- mvrnorm(n = N+d, rep(0, p), sigma)
  X <- matrix(0,N+d,p)
  if(is.null(initX)){
    for(i in 1:d){
      X[i,] <- epi[i,]
    }
  }else{
    for(i in 1:d){
      X[i,] <- initX[(i-1)*p+1,i*p]
    }
  }
  
  for(i in (d+1):(N+d)){
    for(j in 1:d){
      X[i,] <- A[,((j-1)*p+1):(j*p)]%*%X[i-j,] + X[i,]
    }
    X[i,] <- X[i,]+epi[i,]
  }
  X[(d+1):(N+d),]
}



# the following function is only for plot 
Mat_img <- function(beta, len){
  sign(sign(beta - len)*sign(beta+len)+1)
}

image_channel <- function(mat, zlim= c(0,1), ax = NULL) { 
  mat <- t(mat)[,nrow(mat):1]
  image(mat, zlim = zlim, axes = FALSE,col = gray(seq(1, 0, length = 256)))
  box()
  if(!is.null(ax)){
    axis(1, at = seq(0, 1, length = nrow(mat)), labels = rownames(mat))
    axis(2, at = seq(0, 1, length = ncol(mat)), labels = colnames(mat))  
  }
}

image_segment <- function(mat){ 
  mat <- t(mat)[,nrow(mat):1]
  image(mat, axes = FALSE,col = heat.colors(12))
  axis(1, at = seq(0, 1, length = nrow(mat)), labels = rownames(mat))
  axis(2, at = seq(0, 1, length = ncol(mat)), labels = colnames(mat))
  box() 
}

image_segment <- function(mat){
  n <- nrow(mat)
  
  xmin = 0; xmax = 1
  #Generate the palette for the matrix and the legend.  Generate labels for the legend
  palmat <- color.scale(mat, c(1, 0.4), c(1, 0.4), c(0.96, 1))
  palleg <- color.gradient(c(1, 0.4), c(1, 0.4), c(0.96, 1), nslices=100)
  lableg <- c(formatC(xmin, format="f", digits=2), formatC(1*(xmax-xmin)/4, format="f", digits=2), formatC(2*(xmax-xmin)/4, format="f", digits=2), formatC(3*(xmax-xmin)/4, format="f", digits=2), formatC(xmax, format="f", digits=2))
  
  #Set up the plot area and plot the matrix
  par(mar=c(3, 3, 3, 3))
  color2D.matplot(mat, cellcolors=palmat, show.values=2, vcol=rgb(0,0,0), axes=FALSE, vcex=0.7)
  axis(1, at=seq(1, n, 1)-0.5, labels=seq(1, n, 1), tck=-0.01, padj=-1)
  
  #In the axis() statement below, note that the labels are decreasing.  This is because
  #the above color2D.matplot() statement has "axes=FALSE" and a normal axis()
  #statement was used.
  axis(2, at=seq(1, n, 1)-0.5, labels=seq(n, 1, -1), tck=-0.01, padj=0.7)
  
  #Plot the legend
  pardat <- par()
  #color.legend(pardat$usr[2]+0.5, 0, pardat$usr[2]+1, pardat$usr[2], paste(" ", lableg, sep=""), palleg, align="rb", gradient="y", cex=0.7)
}

get_imgmatrix <- function(Mat,d){
  p <- ncol(Mat)
  output<- matrix(0,p,p)
  for(i in 1:d){
    output <- output + Mat[((i-1)*p+1):(i*p),1:p]
  }
  sign(output)
}

find_real <- function(i){
  if(i == 1){
    return (c(2,4,6));
  }else if(i == 2){
    return (c(3,5,7));
  }else if(i == 3){
    return (c(1,4,6));
  }else if(i == 4){
    return (c(2,5,7));
  }else if(i == 5){
    return (c(1,3,6));
  }else if(i == 6){
    return (c(2,4,7));
  }else if(i == 7){
    return (c(1,3,5));
  }else{
    return (c(2,4,6));
  }
}

plot_hamming <- function(ls_mat){
  n <- length(ls_mat)
  p <- nrow(ls_mat[[1]])
  Mat_distance <- matrix(0,n,n)
  for(i in 1:n){
    for(j in 1:n){
      Mat_distance[i,j] <- Hamming_distance(ls_mat[[i]], ls_mat[[j]])/(p^2)
    }
  }
  Mat_distance
}

Hamming_distance <- function(Mat1, Mat2){
  return (sum(abs(Mat1-Mat2)));
}
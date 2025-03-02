

library(magick)
library(plot.matrix)

source("func.R")

#---------------------------------------------------------------------------------------------------------------------
#process data
resize <- function(x1,y1,x2,y2,mat){
  i1 <- x2/x1; i2 <- y2/y1
  mat_new <- matrix(0,x1,y1)
  for(i in 1:x1){
    for(j in 1:y1){
      mat_new[i,j] <- sum(mat[(1+(i-1)*i1):(i*i1),(1+(j-1)*i2):(j*i2)])/(i1*i2)
    }
  }
  return(mat_new)
}


frame <- c()
path <- 'JPEGS/Meet_WalkTogether'
for(i in 2000:2826){
  name <- paste(as.character(i),sep="",".jpg")
  filename <- paste(path, sep = "", name)
  img_temp <- image_read(filename)
  img_gray <-image_convert(img_temp, type = 'grayscale')
  mat0 <- image_data(img_gray, 'rgba')
  mat <- as.integer(mat0)
  mat_new <- as.matrix(mat[1:288,1:384,1],288,384)
  mat_new <- resize(24,32,288,384, mat_new)/256
  mat_vec <- c(mat_new)
  frame <- rbind(frame, mat_vec)
}


#---------------------------------------------------------------------------------------------------------------------
#save and read processed data
# write.matrix(frame,"walk.txt",sep = "\t")
# frame <- as.matrix(read.table("E:/material/matrix completion/walk.txt",sep = "\t", head = FALSE))


#---------------------------------------------------------------------------------------------------------------------
#first step:
#estimate L and S_1,...,S_K
group <- c(1,116,174,232,290,348,406, 464,522,580,638,696)#starting point of each segment
K <- 12 #number of segment
p<- 768; N <- 826
X <- frame[1:826,]
Y <- frame[2:827,]
tau1 <- 0.01; tau2 <- 0.1
lambda <- tau1*sqrt(60*log(768))
mu <- tau2*sqrt(826*768)
niter <- 100
out <- fista.LpS(X, Y, lambda, mu, group, K, niter, backtracking = TRUE)
lr.comp <- out$lr.comp #low rank matrix estimator



#-------------------------------------------------------------------------------------------------------------
#second step: transfer learning for sparse matrix
ls_S_tran <- list();ls_S_selectA <- list() ;ls_S_lasso <- list(); tau <- 0.01
group <- c(group, N+1); ls_select_index <- matrix(0,K,K)
for(i in 1:K){
  n1 <- group[i+1]-group[i]; n2 <- N - n1
  lambda1 = sqrt(tau*log(p)/n1); lambda2 = sqrt(tau*log(p)/N)
  
  target_index <- group[i]:(group[i+1]-1)
  X0 <- X[target_index,]; XA <- X[-target_index,]
  Y0 <- Y[target_index,]; YA <- Y[-target_index,]
  X0_train <- X0[1:floor(2/3*n1),]; X0_test <- X0[(1+floor(2/3*n1)) : n1,]
  Y0_train <- Y0[1:floor(2/3*n1),]; Y0_test <- Y0[(1+floor(2/3*n1)) : n1,]
  X_temp <- rbind(X0_train, XA); Y_temp <- rbind(Y0_train, YA)
  
  #lasso-------------------------------------------------------------------------------------------------------
  ls_S_lasso[[i]] <- lasso(X0_train, Y0_train - X0_train%*%lr.comp, p ,1, lambda1) 
  
  #transfer learning with all data-----------------------------------------------------------------------------
  ls_S_tran[[i]] <- tran_learn(X_temp, Y_temp - X_temp%*%lr.comp, 0.1*lambda1,lambda2, p, 1, floor(2/3*n1), n2)  
  
  #transfer learning with select info set---------------------------------------------------------------------------
  #step 1
  n <- floor(2/3*n1)
  sp = floor(n/2)
  lambda1_test = sqrt(tau*log(p)/sp)
  beta_la_1 = lasso(X0_train[1:sp,], Y0_train[1:sp,] - X0_train[1:sp,]%*%lr.comp, p ,d, lambda1_test)
  beta_la_2 = lasso(X0_train[(sp+1):n,], Y0_train[(sp+1):n,] - X0_train[(sp+1):n,]%*%lr.comp, p ,d, lambda1_test)
  select_set <- c()
  ls_R_k <- c()
  aux_group <- (1:K)[-i]
  for(l in aux_group){
    # X_auxiliary <- XA[((l-1)*(N2-tran)+1):(l*(N2-tran)),]
    X_temp <- rbind(X0_train[1:sp,], X[group[l]:(group[l+1]-1),])
    Y_temp <- rbind(Y0_train[1:sp,], Y[group[l]:(group[l+1]-1),])
    
    lambda2_test = sqrt(tau*log(p)/(sp+(group[l+1]-(group[l]))))
    beta_tran_temp = lasso(X_temp, Y_temp - X_temp%*%lr.comp, p, 1, lambda2_test)
    
    X_test <- X0_train[(sp+1):n,]
    Y_test <- Y0_train[(sp+1):n,]
    
    R_k = sum((Y_test - X_test%*%beta_tran_temp - X_test%*%lr.comp)^2) - sum((Y_test - X_test%*%beta_la_1 - X_test%*%lr.comp)^2) 
    #print(R_k)
    R_thre <- abs(sum((Y_test - X_test%*%beta_la_2 - X_test%*%lr.comp)^2) - sum((Y_test - X_test%*%beta_la_1 - X_test%*%lr.comp)^2)) 
    ls_R_k <- c(ls_R_k,R_k)
  }
  select_set <- which(ls_R_k < 0.1*R_thre)
  ls_select_index[i,aux_group[select_set]] <- 1 
  select_index <- 1:floor(2/3*n1) + group[i] - 1; N <- n1
  for(l in aux_group[select_set]){
    select_index <- c(select_index, group[l]:(group[l+1]-1) )
    N <- N + group[l+1] - group[l]
  }
  
  lambda2 = sqrt(tau*log(p)/(N))
  
  #step2: 
  ls_S_selectA[[i]]  = tran_learn(X[select_index,], Y[select_index,] - X[select_index,]%*%lr.comp , 0.1*lambda1, lambda2, p, 1, floor(2/3*n1), N)
}


#-----------------------------------------------------------------------------------------------------------------------------
#output -- MSE and sd
res_tran <- c(); res_lasso <- c(); res_selectA <- c()
sd_tran <- c(); sd_lasso <- c(); sd_selectA <- c()

for(i in 1:K){
  n1 <- group[i+1]-group[i]; n2 <- N - n1
  lambda1 = sqrt(tau*log(p)/n1); lambda2 = sqrt(tau*log(p)/N)
  
  target_index <- group[i]:(group[i+1]-1)
  X0 <- X[target_index,]; XA <- X[-target_index,]
  Y0 <- Y[target_index,]; YA <- Y[-target_index,]
  X0_train <- X0[1:floor(2/3*n1),]; X0_test <- X0[(1+floor(2/3*n1)) : n1,]
  Y0_train <- Y0[1:floor(2/3*n1),]; Y0_test <- Y0[(1+floor(2/3*n1)) : n1,]
  
  res_tran <- c(res_tran, sum((Y0_test - X0_test%*%lr.comp - X0_test%*%ls_S_tran[[i]])^2))
  res_lasso <- c(res_lasso, sum((Y0_test - X0_test%*%lr.comp - X0_test%*%ls_S_lasso[[i]])^2))
  res_selectA <- c(res_selectA, sum((Y0_test - X0_test%*%lr.comp - X0_test%*%ls_S_selectA[[i]])^2))
  
  sd_tran <- c(sd_tran, sd((Y0_test - X0_test%*%lr.comp - X0_test%*%ls_S_tran[[i]])))
  sd_lasso <- c(sd_lasso, sd((Y0_test - X0_test%*%lr.comp - X0_test%*%ls_S_lasso[[i]])))
  sd_selectA <- c(sd_selectA, sd((Y0_test - X0_test%*%lr.comp - X0_test%*%ls_S_selectA[[i]])))
}



#----------------------------------------------------------------------------------------------------------------------------------
#output -- sparse matrix estimator: S_1,...,S_K
#lasso -- ls_S_lasso;  transfer-lasso -- ls_S_tran;  
ls_MAT<- list()
for(i in 1:K){
  #S_temp <- ls_S_lasso[[i]];
  S_temp <- ls_S_selectA[[i]]
  ni <- group[i+1] - group[i]
  n_sp <- floor((group[i+1] - group[i])*2/3)
  target_index <- group[i]:(group[i+1]-1)
  X0 <- X[target_index,]; XA <- X[-target_index,]
  Y0 <- Y[target_index,]; YA <- Y[-target_index,]
  X0_train <- X0[1:n_sp,]; X0_test <- X0[(1+n_sp) : ni,]
  Y0_train <- Y0[1:n_sp,]; Y0_test <- Y0[(1+n_sp) : ni,]
  Sigma_e_tran = t(Y0_train - X0_train%*%lr.comp - X0_train%*%S_temp)%*%(Y0_train - X0_train%*%lr.comp - X0_train%*%S_temp)/n_sp
  fit_tran_debais = lasso_debais(X0_train, Y0_train - X0_train%*%out$lr.comp, Sigma_e_tran, alpha=2 , S_temp, p, 1, mu =2)
  matimg_temp <- Mat_img(S_temp, qnorm(0.975)*sqrt(fit_tran_debais$V/n_sp))
  
  index <- which(matimg_temp > 0)
  MAT_index <- matrix(0, 768, 768)
  MAT_index[index] <- 1
  MAT <- matrix(apply(MAT_index,2,sum),24,32)
  ls_MAT[[i]] <- MAT
}


#image of sparse matrix
MAX <- 0
for(i in 1:K){
  name <- paste("segment",sep="",as.character(i))
  filename <- paste(name, sep="",".txt")
  #MAT <- as.matrix(read.table(filename,sep = "\t",head = FALSE))
  MAX <- max(c(MAX, max(ls_MAT[[i]])))
  #MAx <- max(c(MAX, max(MAT)))
}

for(i in 1:K){
  name <- paste("segment",sep="",as.character(i))
  filename <- paste(name,sep="",".txt")
  #MAT <- as.matrix(read.table(filename,sep = "\t",head = FALSE))
  MAT <- ls_MAT[[i]]/MAX
  image_channel(MAT)
}


# This file contains main functions for the sparse ISLET.
library(rTensor)
library(matrixcalc)

#' Sparse ISLET
#'
#' This is the main function to do ISLET in the sparse tensor regression case with known rank r, the output is the estimation error and the running time
#' @param X input tensor, n*p1*p2*p3
#' @param y input response y
#' @param r input rank like r1,r2,r3
#' @param p input dimension like p1, p2, p3
#' @param sparse_mode Indicator of the sparse mode of the coefficient tensor, i.e. c(TRUE, TRUE, TRUE)
#' @return Estimate of coefficient tensor A.
#' @examples 
#' set.seed(2018)
#' library(gglasso)
#' # EXAMPLE: Sparse tensor regression
#' n = c(3000); p = c(20); r = c(3); sig = c(5); s = c(12)
#' p1 = p2 = p3 = p
#' r1 = r2 = r3 = r
#' s1 = s2 = s3 = s
#' # Generate data
#' S.array = array(rnorm(r1*r2*r3), dim=c(r1,r2,r3)) # Core tensor
#' S = as.tensor(S.array)
#' E1 = matrix(rnorm(p1*r1), nrow = p1, ncol=r1)
#' E2 = matrix(rnorm(p2*r2), nrow = p2, ncol=r2)
#' E3 = matrix(rnorm(p3*r3), nrow = p3, ncol=r3)
#' # Add sparsity structure
#' E1[(s1+1):p1,] <- 0
#' E2[(s2+1):p2,] <- 0
#' E3[(s3+1):p3,] <- 0
#' sparse_mode = c(TRUE, TRUE, TRUE);
#' X = as.tensor(array(rnorm(n*p1*p2*p3), dim=c(n,p1,p2,p3)))
#' A = ttm(ttm(ttm(S, E1, 1), E2, 2), E3, 3) # Parameter tensor
#' eps = sig * rnorm(n) # Error
#' y = regression.product(X, A) + eps # response
#' snr <- sd(y)/sd(eps)
#' sparse_sketch_result <- sparse_sketch(X, y, sparse_mode, r1, r2, r3, p1, p2, p3)
#' fnorm(sparse_sketch_result$hatA - A)/fnorm(A) # errorss
#' sparse_sketch_result$running_time
#' @referencess
#' \insertRef{zhang2019islet}{ISLET}
sparse_ISLET <- function(X, y, sparse_mode, r1, r2, r3, p1, p2, p3){
  W = as.tensor(array( t(t(y)%*%k_unfold(X,1)@data), dim=c(p1, p2, p3)))/n
  # sigma_est <- get.sigma(W)
  # STATSVD finds the probing direction
  
  # ESTIMATE SIGMA By Tuning
  sigma_try_list <- median(abs(W@data))/qnorm(0.75) * 0.7^(0:5)
  for (ii in 1:length(sigma_try_list)){
    sigma <- sigma_try_list[ii]
    STATSVD_result <- tryCatch( {STATSVD(Y = W, r = c(r1, r2, r3), sigma = sigma, sparse_mode = sparse_mode)}, error = function(e){
    myerror = 1
    result = list(myerror = myerror)
    return(result)
    })
    if (STATSVD_result$myerror == 0){
      hatU1 = STATSVD_result$U_t[[1]]; hatU2 = STATSVD_result$U_t[[2]]; hatU3 = STATSVD_result$U_t[[3]];
      break
    }
  }
  
  ## USE meadian(abs(tilda(A)))/z_0.75 to estimate sigma
  t <- proc.time()
  Z <- ttm(ttm(ttm(W, t(hatU1), 1), t(hatU2), 2), t(hatU3), 3)
  hatV1 = svd(k_unfold(Z, 1)@data)$v
  hatV2 = svd(k_unfold(Z, 2)@data)$v 
  hatV3 = svd(k_unfold(Z, 3)@data)$v

  # to save computation
  myU1 <- ttm(X, t(hatU1), 2)
  myU2 <- ttm(X, t(hatU2), 3)
  rm(X)
  sparse_sketch_result <- sparse_sketch_method(myU1, myU2, hatU1, hatU2, hatU3, hatV1, hatV2, hatV3, r1, r2, r3, p1, p2, p3, y, sparse_mode)
  running_time <- (proc.time() - t)[3]
  result <- list(hatA = sparse_sketch_result$hat_A, running_time = running_time)
  return(result)
}



sparse_sketch_method <- function(myU1, myU2, hatU1, hatU2, hatU3, hatV1, hatV2, hatV3, r1, r2, r3, p1, p2, p3, y, sparse_mode){
  ### INPUT ###
  # myU1, myU2: used to save computation
  # hatU1, hatU2, hatU3: probing direction
  # hatV1, hatV2, hatV3: second sketching direction
  # r1, r2, r3: input rank
  # p1, p2, p3: input dimension
  # y: response
  # sparse_mode: sparse mode of the coefficient tensor
  ### OUTPUT ###
  # estimate of coefficient tensor.

  U2U3 <- ttm(myU2, t(hatU3),4)
  U1U2 <- ttm(myU1, t(hatU2), 3)
  U1U3 <- ttm(myU1, t(hatU3), 4)
  
  # least square to find the Cross estimates
  # Body estimates: Least square
  tildeX_B <- as.matrix(k_unfold(ttm(U2U3, t(hatU1), 2), 1)@data)
  hatB = as.tensor(array(solve((t(tildeX_B)%*%tildeX_B), (t(tildeX_B) %*% y)), c(r1, r2, r3)))

  # Arm estimate: group lasso
  tildeX_E1_intermediate <- ttm(tpartunfold(U2U3, c(3,4)), t(hatV1), 3)
  tildeX_E2_intermediate = ttm(tpartunfold(U1U3, c(2,4)), t(hatV2), 3)
  tildeX_E3_intermediate = ttm(tpartunfold(U1U2, c(2,3)), t(hatV3), 3)
  tildeX_E1_intermediate = k_unfold(tildeX_E1_intermediate, 1)@data
  tildeX_E2_intermediate = k_unfold(tildeX_E2_intermediate, 1)@data
  tildeX_E3_intermediate = k_unfold(tildeX_E3_intermediate, 1)@data

  tildeX_E1 <- array(0, dim = c(n, p1*r1))
  tildeX_E2 <- array(0, dim = c(n, p2*r2))
  tildeX_E3 <- array(0, dim = c(n, p3*r3))

  for (zz in 1:p1){
    selec_index <- seq(from = zz, by = p1, length.out = r1)
    tildeX_E1[,((zz-1)*r1 + 1):(zz*r1)] <- tildeX_E1_intermediate[, selec_index]
  }

  for (zz in 1:p2){
    selec_index <- seq(from = zz, by = p2, length.out = r2)
    tildeX_E2[,((zz-1)*r2 + 1):(zz*r2)] <- tildeX_E2_intermediate[, selec_index]
  }

  for (zz in 1:p3){
    selec_index <- seq(from = zz, by = p3, length.out = r3)
    tildeX_E3[,((zz-1)*r3 + 1):(zz*r3)] <- tildeX_E3_intermediate[, selec_index]
  }
  if (sparse_mode[1]){
    index1 <- rep(1:p1, each = r1)
    E1_data <- data.frame(y = y, tildeX_E1)
    cv1 <- cv.gglasso(tildeX_E1, y, group = index1, pred.loss = "L2", loss = "ls")
    fit1 <- gglasso(tildeX_E1, y, group = index1, loss = "ls", lambda = cv1$lambda.min, intercept = F)
    hat_E1_inter <- fit1$beta
    hat_E1 <- matrix(hat_E1_inter[,1], nrow = p1, ncol = r1, byrow = T)
    } else{
      hat_E1_inter <- solve( (t(tildeX_E1) %*% tildeX_E1 ), (t(tildeX_E1) %*% y)) 
      hat_E1 <- matrix(hat_E1_inter, nrow = p1, ncol = r1, byrow = T)
    }
  
  if (sparse_mode[2]){
      index2 <- rep(1:p2, each = r2)
    E2_data <- data.frame(y = y, tildeX_E2)
    cv2 <- cv.gglasso(tildeX_E2, y, group = index2, pred.loss = "L2", loss = "ls")
    fit2 <- gglasso(tildeX_E2, y, group = index2, loss = "ls", lambda = cv2$lambda.min, intercept = F)
    hat_E2_inter <- fit2$beta
    hat_E2 <- matrix(hat_E2_inter[,1], nrow = p2, ncol = r2, byrow = T)
    } else{
      hat_E2_inter <- solve( (t(tildeX_E2) %*% tildeX_E2 ), (t(tildeX_E2) %*% y)) 
      hat_E2 <- matrix(hat_E2_inter, nrow = p2, ncol = r2, byrow = T)
    }

  if (sparse_mode[3]){
    index3 <- rep(1:p3, each = r3)
    E3_data <- data.frame(y = y, tildeX_E3)
    cv3 <- cv.gglasso(tildeX_E3, y, group = index3, pred.loss = "L2", loss = "ls")
    fit3 <- gglasso(tildeX_E3, y, group = index3, loss = "ls", lambda = cv3$lambda.min, intercept = F)
    hat_E3_inter <- fit3$beta
    hat_E3 <- matrix(hat_E3_inter[,1], nrow = p3, ncol = r3, byrow = T)
  } else{
    hat_E3_inter <- solve( (t(tildeX_E3) %*% tildeX_E3 ), (t(tildeX_E3) %*% y))
    hat_E3 <- matrix(hat_E3_inter, nrow = p3, ncol = r3, byrow = T)
  }
  # tcross.result = tcross(hatB, list((hat_E1 %*% solve(t(hatU1) %*% hat_E1)), (hat_E2 %*% solve(t(hatU2) %*% hat_E2)), (hat_E3 %*% solve(t(hatU3) %*% hat_E3))))
  # hat_A <- tcross.result[[1]]
  hat_A <- ttl(hatB, list((hat_E1 %*% solve(t(hatU1) %*% hat_E1)), (hat_E2 %*% solve(t(hatU2) %*% hat_E2)), (hat_E3 %*% solve(t(hatU3) %*% hat_E3))), c(1,2,3) )
  result <- list(hat_A = hat_A)
  return(result)
}

#' Estimate rank in the sparse ISLET case with unknown rank
#'
#' This function is used to select rank r for the sparse tensor regression case when the rank is unknonw
#' @param y input response y
#' @param r input conservative esimate rank
#' @param p input dimension
#' @param thres a thresholding in the autoselection of rank r. Suggestion value 1-2.
#' @return Estimated rank
#' @examples
#' set.seed(2018)
#' n = c(3000); p = c(20); r = c(3); sig = c(5); s = c(12);
#' thres = c(2) # thresholding is greater than 1, and not too big. Suggestion is 1-2.
#' p1 = p2 = p3 = p
#' r1 = r2 = r3 = r
#' s1 = s2 = s3 = s
#' # Generate data
#' S.array = array(rnorm(r1*r2*r3), dim=c(r1,r2,r3)) # Core tensor
#' S = as.tensor(S.array)
#' E1 = matrix(rnorm(p1*r1), nrow = p1, ncol=r1)
#' E2 = matrix(rnorm(p2*r2), nrow = p2, ncol=r2)
#' E3 = matrix(rnorm(p3*r3), nrow = p3, ncol=r3)
#' # Add sparsity structure
#' E1[(s1+1):p1,] <- 0
#' E2[(s2+1):p2,] <- 0
#' E3[(s3+1):p3,] <- 0
#' sparse_mode = c(TRUE, TRUE, TRUE);
#' X = as.tensor(array(rnorm(n*p1*p2*p3), dim=c(n,p1,p2,p3)))
#' A = ttm(ttm(ttm(S, E1, 1), E2, 2), E3, 3)
#' eps = sig * rnorm(n) # Error
#' y = regression.product(X, A) + eps # response
#' snr <- sd(y)/sd(eps)
#' r1_init = 10; r2_init = 10; r3_init = 10; # an initial conservative estimator of rank
#' r_select = sparse_select_r(X, y,r1, r2, r3, p1, p2, p3, thres, sparse_mode)
#' sparse_sketch_result <- sparse_sketch(X, y, sparse_mode, r_select[1], r_select[2], r_select[3], p1, p2, p3)
#' fnorm(sparse_sketch_result$hatA - A)/fnorm(A) # errorss
#' sparse_sketch_result$running_time
#' @references
#' \insertRef{zhang2019cross}{ISLET}
#'
#' \insertRef{zhang2019islet}{ISLET}
sparse_select_r <- function(X, y, r1, r2, r3, p1, p2, p3, thres, sparse_mode){
   W = as.tensor(array( t(t(y)%*%k_unfold(X,1)@data), dim=c(p1, p2, p3)))/n
  # sigma_est <- get.sigma(W)
  # STATSVD finds the probing direction
  
  ## USE meadian(abs(tilda(A)))/z_0.75 to estimate sigma
  # ESTIMATE SIGMA By Tuning
  sigma_try_list <- median(abs(W@data))/qnorm(0.75) * 0.7^(0:5)
  for (ii in 1:length(sigma_try_list)){
    sigma <- sigma_try_list[ii]
    STATSVD_result <- tryCatch( {STATSVD(W, c(r1, r2, r3), sigma = sigma)}, error = function(e){
    myerror = 1
    result = list(myerror = myerror)
    return(result)
    })
    if (STATSVD_result$myerror == 0){
      hatU1 = STATSVD_result$U_t[[1]]; hatU2 = STATSVD_result$U_t[[2]]; hatU3 = STATSVD_result$U_t[[3]];
      break
    }
  }

  Z <- ttm(ttm(ttm(W, t(hatU1), 1), t(hatU2), 2), t(hatU3), 3)
  hatV1 = svd(k_unfold(Z, 1)@data)$v
  hatV2 = svd(k_unfold(Z, 2)@data)$v 
  hatV3 = svd(k_unfold(Z, 3)@data)$v
  myU1 <- ttm(X, t(hatU1), 2)
  myU2 <- ttm(X, t(hatU2), 3)

  U2U3 <- ttm(myU2, t(hatU3),4)
  U1U2 <- ttm(myU1, t(hatU2), 3)
  U1U3 <- ttm(myU1, t(hatU3), 4)
  
  hatU1.complete = qr.Q(qr(hatU1),complete=TRUE)
  hatU2.complete = qr.Q(qr(hatU2),complete=TRUE)
  hatU3.complete = qr.Q(qr(hatU3),complete=TRUE)

  # least square to find the Cross estimates
  # Body estimates: Least square
  tildeX_B <- as.matrix(k_unfold(ttm(U2U3, t(hatU1), 2), 1)@data)
  hatB = as.tensor(array(solve((t(tildeX_B)%*%tildeX_B), (t(tildeX_B) %*% y)), c(r1, r2, r3)))

  # Arm estimate: group lasso
  tildeX_E1_intermediate <- ttm(tpartunfold(U2U3, c(3,4)), t(hatV1), 3)
  tildeX_E2_intermediate = ttm(tpartunfold(U1U3, c(2,4)), t(hatV2), 3)
  tildeX_E3_intermediate = ttm(tpartunfold(U1U2, c(2,3)), t(hatV3), 3)
  tildeX_E1_intermediate = k_unfold(tildeX_E1_intermediate, 1)@data
  tildeX_E2_intermediate = k_unfold(tildeX_E2_intermediate, 1)@data
  tildeX_E3_intermediate = k_unfold(tildeX_E3_intermediate, 1)@data

  tildeX_E1 <- array(0, dim = c(n, p1*r1))
  tildeX_E2 <- array(0, dim = c(n, p2*r2))
  tildeX_E3 <- array(0, dim = c(n, p3*r3))

  for (zz in 1:p1){
    selec_index <- seq(from = zz, by = p1, length.out = r1)
    tildeX_E1[,((zz-1)*r1 + 1):(zz*r1)] <- tildeX_E1_intermediate[, selec_index]
  }

  for (zz in 1:p2){
    selec_index <- seq(from = zz, by = p2, length.out = r2)
    tildeX_E2[,((zz-1)*r2 + 1):(zz*r2)] <- tildeX_E2_intermediate[, selec_index]
  }

  for (zz in 1:p3){
    selec_index <- seq(from = zz, by = p3, length.out = r3)
    tildeX_E3[,((zz-1)*r3 + 1):(zz*r3)] <- tildeX_E3_intermediate[, selec_index]
  }
  
  index1 <- rep(1:p1, each = r1)
  E1_data <- data.frame(y = y, tildeX_E1)
  cv1 <- cv.gglasso(tildeX_E1, y, group = index1, pred.loss = "L2", loss = "ls")
  fit1 <- gglasso(tildeX_E1, y, group = index1, loss = "ls", lambda = cv1$lambda.min, intercept = F)
  hat_E1_inter <- fit1$beta
  hat_E1 <- matrix(hat_E1_inter[,1], nrow = p1, ncol = r1, byrow = T)
  
  index2 <- rep(1:p2, each = r2)
  E2_data <- data.frame(y = y, tildeX_E2)
  cv2 <- cv.gglasso(tildeX_E2, y, group = index2, pred.loss = "L2", loss = "ls")
  fit2 <- gglasso(tildeX_E2, y, group = index2, loss = "ls", lambda = cv2$lambda.min, intercept = F)
  hat_E2_inter <- fit2$beta
  hat_E2 <- matrix(hat_E2_inter[,1], nrow = p2, ncol = r2, byrow = T)
  
  index3 <- rep(1:p3, each = r3)
  E3_data <- data.frame(y = y, tildeX_E3)
  cv3 <- cv.gglasso(tildeX_E3, y, group = index3, pred.loss = "L2", loss = "ls")
  fit3 <- gglasso(tildeX_E3, y, group = index3, loss = "ls", lambda = cv3$lambda.min, intercept = F)
  hat_E3_inter <- fit3$beta
  hat_E3 <- matrix(hat_E3_inter[,1], nrow = p3, ncol = r3, byrow = T)
  tcross.result = modi_tcross(hatB, list(t(hatU1.complete) %*% hat_E1,t(hatU2.complete) %*% hat_E2 , t(hatU3.complete) %*% hat_E3), thres)
  hatr <- tcross.result[[2]]
  return(hatr)
}



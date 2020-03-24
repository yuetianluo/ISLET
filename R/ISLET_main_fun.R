# This file contains main function for regular ISLET
library(rTensor)
library(matrixcalc)


sketch_method <- function(myU1, myU2 , hatU1, hatU2, hatU3, hatU1_perp, hatU2_perp, hatU3_perp,hatV1, hatV2, hatV3, r1, r2, r3, p1, p2, p3, y){
  ### INPUT ###
  # myU1, myU2: used to save computation
  # hatU1, hatU2, hatU3: probing direction
  # hatU1_perp, hatU2_perp, hatU3_perp: orthogonal complement of the corresponding part
  # hatV1, hatV2, hatV3: second sketching direction
  # r1, r2, r3: input rank
  # p1, p2, p3: input dimension
  # y: response
  ### OUTPUT ###
  # estimate of core tensor coefficient part and the arm coefficient part

  U2U3 <- ttm(myU2, t(hatU3),4)
  U1U3 <- ttm(myU1, t(hatU3), 4)
  U1U2 <- ttm(myU1, t(hatU2), 3)

  # calculate important covariates
  tildeX_B <- as.matrix(k_unfold(ttm(U2U3, t(hatU1), 2), 1)@data)
  tildeX_D1_intermediate <- ttm(tpartunfold(ttm(U2U3, t(hatU1_perp), 2), c(3,4)), t(hatV1), 3)
  tildeX_D2_intermediate = ttm(tpartunfold(ttm(U1U3, t(hatU2_perp), 3), c(2,4)), t(hatV2), 3)
  tildeX_D3_intermediate = ttm(tpartunfold(ttm(U1U2, t(hatU3_perp), 4), c(2,3)), t(hatV3), 3) 
  tildeX_D1 = k_unfold(tildeX_D1_intermediate, 1)@data
  tildeX_D2 = k_unfold(tildeX_D2_intermediate, 1)@data
  tildeX_D3 = k_unfold(tildeX_D3_intermediate, 1)@data

  # calculate the least square estimates for the combined coefficient
  tildeX = cbind(tildeX_B, tildeX_D1, tildeX_D2, tildeX_D3)
  hat_estimate_combine = solve((t(tildeX)%*%tildeX), (t(tildeX) %*% y))
  hatB_alt = as.tensor(array(hat_estimate_combine[1:(r1*r2*r3)], c(r1, r2, r3)))
  # get hatD1, hatD2, hatD3
  hatD1_armalt = matrix(hat_estimate_combine[(r1*r2*r3+1):(r1*r2*r3+(p1-r1)*r1)], nrow=p1-r1, ncol=r1)
  hatD2_armalt = matrix(hat_estimate_combine[(r1*r2*r3+(p1-r1)*r1+1):(r1*r2*r3+(p1-r1)*r1+(p2-r2)*r2)], nrow=p2-r2, ncol=r2)
  hatD3_armalt = matrix(hat_estimate_combine[(r1*r2*r3+(p1-r1)*r1+(p2-r2)*r2+1):
                                               (r1*r2*r3+(p1-r1)*r1+(p2-r2)*r2+(p3-r3)*r3)], nrow=p3-r3, ncol=r3)
  # get the arm part
  hatD1_alt = hatU1 %*% (k_unfold(hatB_alt, 1)@data) %*% hatV1 + hatU1_perp %*% hatD1_armalt
  hatD2_alt = hatU2 %*% (k_unfold(hatB_alt, 2)@data) %*% hatV2 + hatU2_perp %*% hatD2_armalt
  hatD3_alt = hatU3 %*% (k_unfold(hatB_alt, 3)@data) %*% hatV3 + hatU3_perp %*% hatD3_armalt
  result <- list(hatB_alt = hatB_alt, hatD1_alt = hatD1_alt, hatD2_alt=hatD2_alt, hatD3_alt=hatD3_alt)
  return(result)
}

#' ISLET for regular tensor regression, 3-d case
#'
#' This function do ISLET under the regular low-tucker rank tensor regression case. For detail about the procedure, see Zhang et al 2019.
#' @param X input tensor, n*p1*p2*p3
#' @param y input response y
#' @param r input rank like r1,r2,r3
#' @param p input dimension like p1, p2, p3
#' @return Estimate of coefficient tensor A.
#' @examples
#' set.seed(2018)
#' n = c(1000); p = c(10); r = c(3); sig = c(5);
#' p1 = p2 = p3 = p
#' r1 = r2 = r3 = r
#' # Generate data
#' S.array = array(rnorm(r1*r2*r3), dim=c(r1,r2,r3)) # Core tensor
#' S = as.tensor(S.array)
#' E1 = matrix(rnorm(p1*r1), nrow = p1, ncol=r1)
#' E2 = matrix(rnorm(p2*r2), nrow = p2, ncol=r2)
#' E3 = matrix(rnorm(p3*r3), nrow = p3, ncol=r3)
#' X = as.tensor(array(rnorm(n*p1*p2*p3), dim=c(n,p1,p2,p3)))
#' # Parameter tensor
#' A = ttm(ttm(ttm(S, E1, 1), E2, 2), E3, 3)
#' eps = sig * rnorm(n) # Error
#' y = regression.product(X, A) + eps # response
#' snr <- sd(y)/sd(eps)
#' # Estimator
#' A_sketch <- sketch(X, y, r1, r2, r3, p1, p2, p3)
#' error <- fnorm(A_sketch$hatA - A)/ fnorm(A)
#' @references
#' \insertRef{zhang2019islet}{ISLET}
ISLET <- function(X, y, r1, r2, r3, p1, p2, p3){
  t <- proc.time()
  W = as.tensor(array( t(t(y)%*%k_unfold(X,1)@data), dim=c(p1, p2, p3)))/n
  HOOI_result = tucker(W, c(r1,r2,r3)) # HOOI finds the probing directions
  hatU = HOOI_result$U; hatU1 = hatU[[1]]; hatU2 = hatU[[2]]; hatU3 = hatU[[3]];
 # mean(c(sine.theta(hatU1, U1), sine.theta(hatU2, U2), sine.theta(hatU3, U3)))
  # probing direction
  hatV1 = svd(k_unfold(HOOI_result$Z, 1)@data)$v 
  hatV2 = svd(k_unfold(HOOI_result$Z, 2)@data)$v 
  hatV3 = svd(k_unfold(HOOI_result$Z, 3)@data)$v
  
  hatU1_perp = qr.Q(qr(hatU1),complete=TRUE)[,(r1+1):p1]
  hatU2_perp = qr.Q(qr(hatU2),complete=TRUE)[,(r2+1):p2]
  hatU3_perp = qr.Q(qr(hatU3),complete=TRUE)[,(r3+1):p3]
  hatU1.complete = cbind(hatU1, hatU1_perp)
  hatU2.complete = cbind(hatU2, hatU2_perp)
  hatU3.complete = cbind(hatU3, hatU3_perp)
  # to save the computation
  myU1 <- ttm(X, t(hatU1), 2)
  myU2 <- ttm(X, t(hatU2), 3)

  sketch_result2 <- sketch_method(myU1, myU2 , hatU1, hatU2, hatU3, hatU1_perp, hatU2_perp, hatU3_perp, hatV1, hatV2, hatV3, r1, r2, r3, p1, p2, p3, y)
  hatD1_alt = sketch_result2$hatD1_alt
  hatD2_alt = sketch_result2$hatD2_alt
  hatD3_alt = sketch_result2$hatD3_alt
  hatB_alt <- sketch_result2$hatB_alt
  # final estimate
  hatL1 <- hatD1_alt %*% solve(k_unfold(hatB_alt, 1)@data %*% hatV1)
  hatL2 <- hatD2_alt %*% solve(k_unfold(hatB_alt, 2)@data %*% hatV2)
  hatL3 <- hatD3_alt %*% solve(k_unfold(hatB_alt, 3)@data %*% hatV3)
  hatA_alt = ttl(hatB_alt, list(hatL1, hatL2, hatL3), c(1,2,3))
  running_time <- (proc.time() -t)[3]
  result <- list(hatA = hatA_alt, running_time = running_time)
  return(result)
}


ISLET_initial_est <- function(X, y, r1, r2, r3, p1, p2, p3){
  t <- proc.time()
  W = as.tensor(array( t(t(y)%*%k_unfold(X,1)@data), dim=c(p1, p2, p3)))/n
  HOOI_result = tucker(W, c(r1,r2,r3)) # HOOI finds the probing directions
  tildeA <- HOOI_result$est
  running_time <- (proc.time() -t)[3]
  result <- list(tildeA = tildeA, running_time = running_time)
  return(result)
}


matrix_sketch_method <- function(myU1, myU2, hatU1, hatU2, hatU1_perp, hatU2_perp, r, p1, p2, y){
  ### INPUT ###
  # myU1, myU2: to save the compuation
  # hatU1, hatU2: probing direction 1
  # hatU1_perp, hatU2_perp: orthogonal part of probing direction 1
  # hatV1, hatV2: probing direction 2
  ### OUTPUT ###
  # estimate of core tensor coefficient part and the arm coefficient part

  # calculate the importance covariates
  tildeX_B <- as.matrix(k_unfold(ttm(myU2, t(hatU1), 2), 1)@data)
  tildeX_D1_intermediate <- ttm(myU2, t(hatU1_perp), 2)
  tildeX_D2_intermediate = tpartunfold(ttm(myU1, t(hatU2_perp), 3),c(2))
  tildeX_D1 = k_unfold(tildeX_D1_intermediate, 1)@data
  tildeX_D2 = k_unfold(tildeX_D2_intermediate, 1)@data

  # get combined estimate
  tildeX = cbind(tildeX_B, tildeX_D1, tildeX_D2)
  hat_estimate_combine = solve((t(tildeX)%*%tildeX), (t(tildeX) %*% y))
  hatB_alt = array(hat_estimate_combine[1:(r*r)], c(r, r))
  hatD1_armalt = matrix(hat_estimate_combine[(r*r+1):(r*r+(p1-r)*r)], nrow=p1-r, ncol=r)
  hatD2_armalt = matrix(hat_estimate_combine[(r*r+(p1-r)*r+1):(r*r+(p1-r)*r+(p2-r)*r)], nrow=p2-r, ncol=r)
  # get arm estimate
  hatD1_alt = hatU1 %*% hatB_alt + hatU1_perp %*% hatD1_armalt
  hatD2_alt = hatU2 %*% t(hatB_alt) + hatU2_perp %*% hatD2_armalt
  result <- list(hatB_alt = hatB_alt, hatD1_alt = hatD1_alt, hatD2_alt=hatD2_alt)
  return(result)
}


#' Matrix ISLET, 2-d case
#'
#' @param X input tensor data, n*p1*p2
#' @param y input response y
#' @param r input rank like r
#' @param p input dimension like p1, p2
#' @return Estimate of coefficient tensor A.
#' @examples
#' # EXAMPLE ISLET in Matrix case
#' n = c(5000); p = c(50); r = c(3); sig = c(5);
#' # data generation
#' S = diag(abs(rnorm(r))) # Core matrix
#' E1 = matrix(rnorm(p1*r), nrow = p1, ncol=r)
#' E2 = matrix(rnorm(p2*r), nrow = p2, ncol=r)
#' X = as.tensor(array(rnorm(n*p1*p2), dim=c(n,p1,p2)))
#' A = E1 %*% S %*% t(E2)
#' eps = sig * rnorm(n) # Error
#' y = matrix.regression.product(X, A) + eps # response
#' # Estimation
#' A_sketch <- matrix_sketch(X, y, r, p1, p2)
#' hatA <- A_sketch[[1]]
#' norm(hatA - A, type = "F")/norm(A, type = "F") # error
matrix_ISLET <- function(X, y, r, p1, p2){
  t <- proc.time()
  W = array( t(t(y)%*%k_unfold(X,1)@data), dim=c(p1, p2))/n
  svd_result = svd(W) # SVD finds the probing directions
  hatU1 = svd_result$u[,1:r]; hatU2 = svd_result$v[,1:r];
  # mean(c(sine.theta(hatU1, U1), sine.theta(hatU2, U2)))
  
  hatU1_perp = qr.Q(qr(hatU1),complete=TRUE)[,(r+1):p1]
  hatU2_perp = qr.Q(qr(hatU2),complete=TRUE)[,(r+1):p2]
  
  hatU1.complete = cbind(hatU1, hatU1_perp)
  hatU2.complete = cbind(hatU2, hatU2_perp)
  myU1 <- ttm(X, t(hatU1), 2)
  myU2 <- ttm(X, t(hatU2), 3)
  rm(X)
  matrix_sketch_result <- matrix_sketch_method(myU1, myU2, hatU1, hatU2, hatU1_perp, hatU2_perp, r, p1, p2, y)
  hatD1_alt = matrix_sketch_result$hatD1_alt
  hatD2_alt = matrix_sketch_result$hatD2_alt
  hatB_alt <- matrix_sketch_result$hatB_alt
  # final estimate
  hatL1 <- hatD1_alt %*% solve(hatB_alt)
  hatL2 <- hatD2_alt %*% solve(t(hatB_alt))
  hatA_alt = hatL1 %*% hatB_alt %*% t(hatL2)
  running_time <- (proc.time() -t)[3]
  result <- list(hatA = hatA_alt, running_time = running_time)
  return(result)
}

#' Estimate Rank when rank is unknown- 3-d case
#'
#' This function is used to estimate the rank r when the Tucker rank of tensor is unknown. For detail about the procedure, see reference Zhang et al 2019 and Zhang 2018 CROSS: Efficient low rank tensor completion.s
#' @param X input tensor, n*p1*p2*p3
#' @param y input response y
#' @param r input a conservative rank estimation r1_init,r2_init,r3_inits
#' @param p input dimension like p1, p2, p3
#' @param thres thresholding in the autoselection of rank r, suggestion 1-2.
#' @return Estimate Tucker rank of the parameter tensor
#' @examples
#' # ISLET in the unknown rank setting 
#' n = c(1000); p = c(10); r = c(3); sig = c(5);
#' p1 = p2 = p3 = p
#' r1 = r2 = r3 = r
#' thres = 1.5; # thresholding is greater than 1, and not too big. Suggestion is 1-2.
#' # data generate
#' A = ttm(ttm(ttm(S, E1, 1), E2, 2), E3, 3)
#' eps = sig * rnorm(n) # Error
#' y = regression.product(X, A) + eps # response
#' snr <- sd(y)/sd(eps)
#' # Get a estimate of rank
#' r1_init = r2_init = r3_init = 8
#' r_select <- select_r(X, y, thres, r1_init, r2_init, r3_init, p1, p2, p3)
#' # Do estimation
#' A_sketch <- sketch(X, y, r_select[1], r_select[2], r_select[3], p1, p2, p3)
#' #error
#' fnorm(A_sketch - A)/fnorm(A)
#' @references
#' \insertRef{zhang2019cross}{ISLET}
#'
#' \insertRef{zhang2019islet}{ISLET}
select_r <- function(X, y, thres, r1, r2, r3, p1, p2, p3){
  W = as.tensor(array( t(t(y)%*%k_unfold(X,1)@data), dim=c(p1, p2, p3)))/n
  HOOI_result = tucker(W, c(r1,r2,r3)) # HOOI finds the probing directions
  hatU = HOOI_result$U; hatU1 = hatU[[1]]; hatU2 = hatU[[2]]; hatU3 = hatU[[3]];
  # mean(c(sine.theta(hatU1, U1), sine.theta(hatU2, U2), sine.theta(hatU3, U3)))
  hatV1 = svd(k_unfold(HOOI_result$Z, 1)@data)$v 
  hatV2 = svd(k_unfold(HOOI_result$Z, 2)@data)$v 
  hatV3 = svd(k_unfold(HOOI_result$Z, 3)@data)$v
  
  hatU1_perp = qr.Q(qr(hatU1),complete=TRUE)[,(r1+1):p1]
  hatU2_perp = qr.Q(qr(hatU2),complete=TRUE)[,(r2+1):p2]
  hatU3_perp = qr.Q(qr(hatU3),complete=TRUE)[,(r3+1):p3]
  hatU1.complete = cbind(hatU1, hatU1_perp)
  hatU2.complete = cbind(hatU2, hatU2_perp)
  hatU3.complete = cbind(hatU3, hatU3_perp)
  myU1 <- ttm(X, t(hatU1), 2)
  myU2 <- ttm(X, t(hatU2), 3)
  sketch_result2 <- sketch_method(myU1, myU2 , hatU1, hatU2, hatU3, hatU1_perp, hatU2_perp, hatU3_perp, hatV1, hatV2, hatV3, r1, r2, r3, p1, p2, p3, y)
  hatD1_alt = sketch_result2$hatD1_alt
  hatD2_alt = sketch_result2$hatD2_alt
  hatD3_alt = sketch_result2$hatD3_alt
  hatB_alt <- sketch_result2$hatB_alt
  # final estimate
  hatL1 <- hatD1_alt %*% solve(k_unfold(hatB_alt, 1)@data %*% hatV1)
  hatL2 <- hatD2_alt %*% solve(k_unfold(hatB_alt, 2)@data %*% hatV2)
  hatL3 <- hatD3_alt %*% solve(k_unfold(hatB_alt, 3)@data %*% hatV3)
  #hatA_alt = ttl(hatB_alt, list(hatL1, hatL2, hatL3), c(1,2,3))
  tcross.result = modi_tcross(hatB_alt, list(t(hatU1.complete) %*% hatD1_alt, t(hatU2.complete) %*% hatD2_alt, t(hatU3.complete) %*% hatD3_alt), thres)
  hatr <- tcross.result[[2]]
  return(hatr)
}

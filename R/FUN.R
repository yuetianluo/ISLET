
# Low-rank matrix and tensor regression via sketching
library(rTensor)
library(matrixcalc)
library(MASS)


##################################################
# Additional functions
##################################################
regression.product = function(X, A){
  # X: order-(d+1) tensor, first mode is of dimension n 
  # A: order-d tensor, the dimensions matches all but the first dimensions of X
  y = k_unfold(X, 1)@data %*% as.vector(A@data)
  return(y)
}
matrix.regression.product = function(X, A){
  # X: order-(d+1) tensor, first mode is of dimension n 
  # A: order-d tensor, the dimensions matches all but the first dimensions of X
  y = k_unfold(X, 1)@data %*% as.vector(A)
  return(y)
}
####################################################################
# Sine-theta distance between two singular vector spaces U and hatU
# Default: Frobenius sine-theta distance
####################################################################
sine.theta <- function(U, hatU, q){ # Sine-theta distance between two singular subspaces. 
  # U and hatU should be of the same dimensions
  try(if(missing("q")) q = 2) 
  try(if(nrow(U)!=nrow(hatU) | ncol(U)!=ncol(hatU)) stop("Matrix does not match") )
  
  r = ncol(U)
  v = 1 - (svd(t(U) %*% hatU)$d)^2
  if(is.infinite(q))
    return(max(v))
  else
    return((sum(v^q))^(1/q))
}
##################################################

tperm = function(X, v){
  # permute the modes of X according to v
  return(as.tensor(aperm(as.array(X), v)))
}

tpartunfold = function(X, v){
  # unfold the modes in V for X
  # Example: X: 5*3*4*10 tensor; tpartunfold(X, c(2,3)) yields 5*10*12 tensor
  p = dim(X)
  d = length(p)
  v_hat = c(setdiff(1:d, v), v)
  Y = array(aperm(X@data, v_hat), dim = c(p[setdiff(1:d, v)], prod(p[v])))
  return(as.tensor(Y))
}


########## Cross Tensor Completion ###########
mynorm = function(x, type){if(type=="2"){spectral.norm(x)}else{norm(x,type)}}
tcross <- function(B, Arm, Indices, c_t){
  # Cross: efficient tensor completion
  
  # B: body measurements; tensor class input; required
  # arm: list of d matrices representing arm measurements
  # Indices: list of d sets of indices for locations of observed core tensor, defaut value: the top {1, ..., r} 
  
  try(if(missing("B")) stop("missing argument: B is required as tensor type."))
  try(if(missing("Arm")) stop("invalid input: Arm is required as list of matrices."))
  try(if(missing("c_t")) c_t = 3)
  try(if(class(B) != "Tensor") stop("invalid input: S should be of tensor type."))
  try(if(! is.list(Arm)) stop("invalid input: Arm should be of a list of matrices."))
  m = dim(B)
  d = length(m)
  try(if(d != length(Arm)) stop("invalid input: order of S and arm are incompatible."))
  
  
  g = rep(0, d); p = rep(0, d) #Second step.
  for (t in 1:d){
    try(if(! is.matrix(Arm[[t]])) stop("invalid input: Arm should be of a list of matrices."))
    p[t] = dim(Arm[[t]])[1]; g[t] = dim(Arm[[t]])[2]
    try(if(p[t] < m[t]) stop("invalid input: core tensor and arm matrix's dimentions are incompatible."))
  }
  if(missing("Indices")) {
    Indices = list();
    for (t in 1:d){
      Indices = c(Indices, list(1:m[t]))
    }
  }
  
  hatR = list()  # Initialization for expanding matrix hatR[]
  hatr=rep(0, d)
  for (t in 1:d)
  {
    arm_matrix = Arm[[t]]
    body_matrix = k_unfold(B, t)@data
    U_body = svd(body_matrix)$u
    joint_matrix = arm_matrix[Indices[[t]],]
    V_arm = svd(arm_matrix)$v
    joint_matrix = t(U_body) %*% joint_matrix %*% V_arm
    arm_matrix = arm_matrix %*% V_arm
    
    for (s in min(g[t], m[t]):1){
      if(!is.singular.matrix(as.matrix(joint_matrix[1:s,1:s]))){
        if (mynorm(arm_matrix[, 1:s] %*% ginv(joint_matrix[1:s, 1:s]), "2") <= c_t*sqrt(p[t])/(sqrt(m[t]))){
          hatr[t] = s; break;
        }
      }
    }
    if (hatr[t] >0 ){
      hatRt = arm_matrix[, 1:hatr[t]] %*% ginv(joint_matrix[1:hatr[t], 1:hatr[t]]) %*% t(U_body[,1:hatr[t]])
    } else {
      hatRt = matrix(0, p[t], m[t]);
    }
    hatR = c(hatR, list(hatRt))
  }
  
  hatA = B;
  for (t in 1:d){
    hatA = ttm(hatA, hatR[[t]], t)
  }
  
  return(list(hatA, hatr, hatR))
}
### End of the Algorithm ###



# Modified tcross function
modi_tcross <- function(B, Arm, thres,Indices){
  # Cross: efficient tensor completion
  
  # B: body measurements; tensor class input; required
  # arm: list of d matrices representing arm measurements
  # Indices: list of d sets of indices for locations of observed core tensor, defaut value: the top {1, ..., r} 
  
  try(if(missing("B")) stop("missing argument: B is required as tensor type."))
  try(if(missing("Arm")) stop("invalid input: Arm is required as list of matrices."))
  try(if(class(B) != "Tensor") stop("invalid input: S should be of tensor type."))
  try(if(! is.list(Arm)) stop("invalid input: Arm should be of a list of matrices."))
  m = dim(B)
  d = length(m)
  try(if(d != length(Arm)) stop("invalid input: order of S and arm are incompatible."))
  
  
  g = rep(0, d); p = rep(0, d) #Second step.
  for (t in 1:d){
    try(if(! is.matrix(Arm[[t]])) stop("invalid input: Arm should be of a list of matrices."))
    p[t] = dim(Arm[[t]])[1]; g[t] = dim(Arm[[t]])[2]
    try(if(p[t] < m[t]) stop("invalid input: core tensor and arm matrix's dimentions are incompatible."))
  }
  if(missing("Indices")) {
    Indices = list();
    for (t in 1:d){
      Indices = c(Indices, list(1:m[t]))
    }
  }
  
  hatR = list()  # Initialization for expanding matrix hatR[]
  hatr=rep(0, d)
  for (t in 1:d)
  {
    arm_matrix = Arm[[t]]
    body_matrix = k_unfold(B, t)@data
    U_body = svd(body_matrix)$u
    joint_matrix = arm_matrix[Indices[[t]],]
    V_arm = svd(arm_matrix)$v
    joint_matrix = t(U_body) %*% joint_matrix %*% V_arm
    arm_matrix = arm_matrix %*% V_arm
    
    for (s in min(g[t], m[t]):1){
      if(!is.singular.matrix(as.matrix(joint_matrix[1:s,1:s]))){
        if (mynorm(arm_matrix[, 1:s] %*% ginv(joint_matrix[1:s, 1:s]), "2") <= thres){
          hatr[t] = s; break;
        }
      }
    }
    hatr[hatr == 0] <- 3
    if (hatr[t] >0 ){
      hatRt = arm_matrix[, 1:hatr[t]] %*% ginv(joint_matrix[1:hatr[t], 1:hatr[t]]) %*% t(U_body[,1:hatr[t]])
    } else {
      hatRt = matrix(0, p[t], m[t]);
    }
    hatR = c(hatR, list(hatRt))
  }
  
  hatA = B;
  for (t in 1:d){
    hatA = ttm(hatA, hatR[[t]], t)
  }
  
  return(list(hatA, hatr, hatR))
}



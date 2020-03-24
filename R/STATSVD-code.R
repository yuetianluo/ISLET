# Source code for STATSVD

library(rTensor)
library(MASS)

vector_sum_square <- function(x){
  return(sum(x^2))
}

##############################################
# Estimate sigma with MAD.
##############################################
get.sigma <- function(Y){
  if(missing("Y")) stop("missing argument: Y is required as the tensor type.")
  try(if(class(Y) != "Tensor") stop("invalid input: Y should be of tensor type."))
  X = k_unfold(Y,1)@data
  V = as.vector(X)
  return(1.4826*median(abs(V)))
}

##############################################
# Estimate Tucker rank.
##############################################
get.rank <- function(Y, sigma, sparse_mode){
  if(missing("Y")) stop("missing argument: Y is required as the tensor type.")
  if(missing("sigma")) stop("missing argument: sigma is required as the non-negative double")
  try(if(class(Y) != "Tensor") stop("invalid input: Y should be of tensor type."))
  try(if(sigma <= 0) stop("invalid input: sigma must be positive."))
  p = dim(Y)
  d = length(p)
  r = rep(0,d)
  try(if(missing("sparse_mode")) sparse_mode=rep(TRUE, d))
  p.prod = prod(p)
  p.minus.k = prod(p) / p
  I_0 = list()
  for(i in 1:d){ # Initialization step 1: find significant index set I_k
    Y_i = k_unfold(Y, i)@data
    Y_i_row_norm = apply(Y_i, 1, vector_sum_square)
    #consider sparse model!
    if(isTRUE(sparse_mode[i])){
      I_0 = c(I_0, list((apply(Y_i, 1, vector_sum_square) > sigma^2*(p.minus.k[i]+2*sqrt(p.minus.k[i]*log(p.prod))+2*log(p.prod)))
                        |   (apply(abs(Y_i), 1, max) > 2*sigma*sqrt(log(p.prod)))))
    }
    else
      I_0 = c(I_0, list(rep(TRUE,p[i])))
    # select the significant indices
  }
  tilde_Y = Y
  s = rep(0,d)
  for (i in 1:d){ # Initialization step 2: construct tilde_Y
    tilde_Y = ttm(tilde_Y, diag(I_0[[i]]*1), i)
    s[i] = sum(I_0[[i]]*1)
  }
  
  for(i in 1:d){
    MY = k_unfold(tilde_Y,d)@data
    p.k = dim(MY)[1]
    p.minus.k = dim(MY)[2]
    s.value = svd(MY)$d
    ci = s[i]
    cj = prod(s)/s[i]
    delta = sqrt(ci) + sqrt(cj) + sqrt(2*ci*(1+log(p.k/ci)) + 2*cj*(1+log(p.minus.k/cj))
                                       + 4*log(p.k))
    r[i] = length(which(s.value >= sigma*delta))
    
  }
  return(r)
}


######################################################
# Sparse tensor alternating truncation for SVD.
######################################################
STATSVD <- function(Y, r, sigma, tmax, sparse_mode, vartol){
  
  if(missing("Y")) stop("missing argument: Y is required as the tensor type.")
  if(class(Y) != "Tensor") stop("invalid input: Y should be of tensor type.")
  if(missing("sigma")){
    sigma = get.sigma(Y)
  }
  if(class(sigma) != "numeric") stop("invalid input: sigma should be of numeric type.")
  
  if(missing("r")){
    r = get.rank(Y,sigma)
  }
  if(class(r) != "numeric") stop("invalid input: r should be of numeric type.")
  
  p = dim(Y)
  d = length(p)
  
  if(missing("sparse_mode")) sparse_mode=rep(TRUE, d)
  if(missing("tmax")) tmax = 10
  if(missing("vartol")) vartol = 1e-6
  
  if(is.atomic(r) && length(r)==1){
    r = rep(r, d)
  }
  if(d != length(r)) stop("invalid input: r and the order of Y is incompatible.")
  
  
  p.prod = prod(p)
  p.minus.k = prod(p) / p
  r.minus.k = prod(r) / r # Introduce r_{-k}

  U_t = list(); I_0 = list()
  for(i in 1:d){ # Initialization step 1: find significant index set I_k
    Y_i = k_unfold(Y, i)@data
    Y_i_row_norm = apply(Y_i, 1, vector_sum_square)
    if(isTRUE(sparse_mode[i])){ # sparse mode:
      I_0 = c(I_0, list((apply(Y_i, 1, vector_sum_square) > (sigma^2*(p.minus.k[i]+2*sqrt(p.minus.k[i]*log(p.prod))+2*log(p.prod))))
                        | (apply(abs(Y_i), 1, max) > 2*sigma*sqrt(log(p.prod)))))
    }
    else # non-sparse mode
      I_0 = c(I_0, list(rep(TRUE,p[i])))
  }
  tilde_Y = Y
  for (i in 1:d){ # Initialization step 2: construct tilde_Y
    tilde_Y = ttm(tilde_Y, diag(I_0[[i]]*1), i)
  }
  for (i in 1:d){ # Initialization step 3: find loading U_0k
    if(isTRUE(sparse_mode[i])){ # sparse mode:
      Ui = matrix(0, nrow = p[i], ncol = r[i])
      datai = k_unfold(tilde_Y, i)@data
      selec = (1:p[i])[I_0[[i]]]
      if(length(selec) < r[i]){  # in case that sk < rk
        warning("Rank is larger than support!")
        Ui = (svd(datai)$u[,1:r[i]])
      }
      else{
        Ui[I_0[[i]],] = svd(datai[I_0[[i]],])$u[,1:r[i]]
      }
      U_t = c(U_t, list(t(Ui)))
    }
    else{ # non-sparse mode
      U_t = c(U_t, list(t(svd(k_unfold(tilde_Y, i)@data)$u[,1:r[i]])))
    }
    
  }
  
  
  t = 1; approx = -1;
  svector = 0; # store the accumulative variance
  while(t<tmax){ # Stop criterion: convergence or maximum number of iteration reached
    for(i in 1:d){
      A = ttl(Y, U_t[-i], (1:d)[-i])
      A_matrix = k_unfold(A, i)@data
      if(!isTRUE(sparse_mode[i])){ # non-sparse mode
        svd.result = svd(A_matrix)
        U_t[[i]] = t(svd.result$u[,1:r[i]])
        svector = svd.result$d[1:r[i]]
      }
      else{ # sparse mode: first thresholding:
        A_k_row_norm = apply(A_matrix, 1, vector_sum_square)
        I_k = A_k_row_norm > sigma^2 * (r.minus.k[i] + 2*(sqrt(r.minus.k[i]*log(p.prod))+log(p.prod)))
        selec = (1:p[i])[I_k]
        
        B_matrix = A_matrix[I_k,]
        if(length(selec) < r[i]){ # in case that sk < rk
          warning("Rank is larger than support!")
          hat.V = svd(rbind(B_matrix, matrix(0, r[i]-length(selec))))$v[,1:r[i]]
        }
        else{
          hat.V = svd(B_matrix)$v[,1:r[i]]
        }
        
        # second thresholding:
        bar.A = A_matrix %*% hat.V
        bar.B = matrix(0, nrow(bar.A), ncol(bar.A))
        bar.I_k = apply(bar.A, 1, vector_sum_square) > sigma^2 * (r[i] + 2*(sqrt(r[i]*log(p.prod))+log(p.prod)))
        bar.B[bar.I_k, ] = bar.A[bar.I_k,]
        
        if(length((1:p[i])[bar.I_k]) < r[i]){# in case that sk < rk
          warning("Rank is larger than support!")
          svd.result = svd(bar.B)
          svector = svd.result$d
          U_t[[i]] = t(svd.result$u)
        }
        else{
          This.U = matrix(0, nrow(bar.A), r[i])
          C = bar.A[bar.I_k,]
          svd.result = svd(C)
          This.U[bar.I_k,] = svd.result$u[,1:r[i]]
          U_t[[i]] = t(This.U) 
          svector = svd.result$d[1:r[i]]
        }
        
        
      }
    }
    if (abs(sum(svector^2) - approx) > vartol & t<tmax){
      #print(t)
      t = t+1
      approx = sum(svector^2)
    }
    else {
      break
    }
  }
  #print(bar.I_k)
  for(i in 1:d){
    U_t[[i]] = t(U_t[[i]])
  }
  result <- list(U_t = U_t, myerror = 0)
  return(result)
}
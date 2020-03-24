# EXAMPLE1: regular tensor regressison

library(rTensor)
library(matrixcalc)
set.seed(2018)
n = c(2000); p = c(10); r = c(3); sig = c(3);

p1 = p2 = p3 = p
r1 = r2 = r3 = r
# Generate data
S.array = array(rnorm(r1*r2*r3), dim=c(r1,r2,r3)) # Core tensor
S = as.tensor(S.array)
E1 = matrix(rnorm(p1*r1), nrow = p1, ncol=r1)
E2 = matrix(rnorm(p2*r2), nrow = p2, ncol=r2)
E3 = matrix(rnorm(p3*r3), nrow = p3, ncol=r3)
X = as.tensor(array(rnorm(n*p1*p2*p3), dim=c(n,p1,p2,p3)))
# Parameter tensor
A = ttm(ttm(ttm(S, E1, 1), E2, 2), E3, 3)
eps = sig * rnorm(n) # Error
y = regression.product(X, A) + eps # response
snr <- sd(y)/sd(eps)

# Estimator
A_sketch <- ISLET(X, y, r1, r2, r3, p1, p2, p3)
A_sketch_initial <- ISLET_initial_est(X, y, r1, r2, r3, p1, p2, p3)
fnorm(A_sketch$hatA - A)/ fnorm(A) # error
A_sketch$running_time # runnning time
fnorm(A_sketch_initial$tildeA - A)/ fnorm(A) # error
A_sketch_initial$running_time # runnning time


# EXAMPLE2: Regular tensor regression when rank is unknown
set.seed(2010)
n = c(1000); p = c(10); r = c(3); sig = c(3);
p1 = p2 = p3 = p
r1 = r2 = r3 = r
thres = 1.5; # thresholding is greater than 1, and not too big. Suggestion is 1-3.
# data generate
A = ttm(ttm(ttm(S, E1, 1), E2, 2), E3, 3)
eps = sig * rnorm(n) # Error
y = regression.product(X, A) + eps # response
snr <- sd(y)/sd(eps)

# Get a estimate of ranks
r1_init = r2_init = r3_init = 8
r_select <- select_r(X, y, thres, r1_init, r2_init, r3_init, p1, p2, p3)
r_selectss
# Do estimation
A_sketch <- ISLET(X, y, r_select[1], r_select[2], r_select[3], p1, p2, p3)
#error
fnorm(A_sketch$hatA - A)/fnorm(A)
A_sketch$running_time

set.seed(2018)
library(gglasso)
# EXAMPLE3: Sparse tensor regression
n = c(3000); p = c(20); r = c(3); sig = c(5); s = c(12)
p1 = p2 = p3 = p
r1 = r2 = r3 = r
s1 = s2 = s3 = s
# Generate data
S.array = array(rnorm(r1*r2*r3), dim=c(r1,r2,r3)) # Core tensor
S = as.tensor(S.array)
E1 = matrix(rnorm(p1*r1), nrow = p1, ncol=r1)
E2 = matrix(rnorm(p2*r2), nrow = p2, ncol=r2)
E3 = matrix(rnorm(p3*r3), nrow = p3, ncol=r3)
# Add sparsity structure
E1[(s1+1):p1,] <- 0
E2[(s2+1):p2,] <- 0
E3[(s3+1):p3,] <- 0
sparse_mode = c(TRUE, TRUE, TRUE);
X = as.tensor(array(rnorm(n*p1*p2*p3), dim=c(n,p1,p2,p3)))

A = ttm(ttm(ttm(S, E1, 1), E2, 2), E3, 3) # Parameter tensor
eps = sig * rnorm(n) # Error
y = regression.product(X, A) + eps # response
snr <- sd(y)/sd(eps)

sparse_sketch_result <- sparse_ISLET(X, y, sparse_mode, r1, r2, r3, p1, p2, p3)
fnorm(sparse_sketch_result$hatA - A)/fnorm(A) # errorss
sparse_sketch_result$running_time


# EXAMPLE4: Sparse ISLET when the rank is unknwon
set.seed(2018)
n = c(3000); p = c(20); r = c(3); sig = c(5); s = c(12);
thres = c(2) # thresholding is greater than 1, and not too big. Suggestion is 1-3.
p1 = p2 = p3 = p
r1 = r2 = r3 = r
s1 = s2 = s3 = s
# Generate data
S.array = array(rnorm(r1*r2*r3), dim=c(r1,r2,r3)) # Core tensor
S = as.tensor(S.array)
E1 = matrix(rnorm(p1*r1), nrow = p1, ncol=r1)
E2 = matrix(rnorm(p2*r2), nrow = p2, ncol=r2)
E3 = matrix(rnorm(p3*r3), nrow = p3, ncol=r3)
# Add sparsity structure
E1[(s1+1):p1,] <- 0
E2[(s2+1):p2,] <- 0
E3[(s3+1):p3,] <- 0
sparse_mode = c(TRUE, TRUE, TRUE);
X = as.tensor(array(rnorm(n*p1*p2*p3), dim=c(n,p1,p2,p3)))

A = ttm(ttm(ttm(S, E1, 1), E2, 2), E3, 3)
eps = sig * rnorm(n) # Error
y = regression.product(X, A) + eps # response
snr <- sd(y)/sd(eps)

r1_init = 10; r2_init = 10; r3_init = 10; # an initial conservative estimator of rank
r_select = sparse_select_r(X, y,r1, r2, r3, p1, p2, p3, thres, sparse_mode)
sparse_sketch_result <- sparse_ISLET(X, y, sparse_mode, r_select[1], r_select[2], r_select[3], p1, p2, p3)
fnorm(sparse_sketch_result$hatA - A)/fnorm(A) # errorss
sparse_sketch_result$running_time

# EXAMPLE5 ISLET in Matrix case
n = c(5000); p = c(50); r = c(3); sig = c(5);
# data generation
S = diag(abs(rnorm(r))) # Core matrix
E1 = matrix(rnorm(p1*r), nrow = p1, ncol=r)
E2 = matrix(rnorm(p2*r), nrow = p2, ncol=r)
X = as.tensor(array(rnorm(n*p1*p2), dim=c(n,p1,p2)))
A = E1 %*% S %*% t(E2)
eps = sig * rnorm(n) # Error
y = matrix.regression.product(X, A) + eps # response

# Estimation
A_sketch <- matrix_ISLET(X, y, r, p1, p2)
hatA <- A_sketch[[1]]
norm(hatA - A, type = "F")/norm(A, type = "F") # error

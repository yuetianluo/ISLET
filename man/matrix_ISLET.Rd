% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ISLET_main_fun.R
\name{matrix_ISLET}
\alias{matrix_ISLET}
\title{Matrix ISLET, 2-d case}
\usage{
matrix_ISLET(X, y, r, p1, p2)
}
\arguments{
\item{X}{input tensor data, n*p1*p2}

\item{y}{input response y}

\item{r}{input rank like r}

\item{p}{input dimension like p1, p2}
}
\value{
Estimate of coefficient tensor A.
}
\description{
Matrix ISLET, 2-d case
}
\examples{
# EXAMPLE ISLET in Matrix case
n = c(5000); p = c(50); r = c(3); sig = c(5);
# data generation
S = diag(abs(rnorm(r))) # Core matrix
E1 = matrix(rnorm(p1*r), nrow = p1, ncol=r)
E2 = matrix(rnorm(p2*r), nrow = p2, ncol=r)
X = as.tensor(array(rnorm(n*p1*p2), dim=c(n,p1,p2)))
A = E1 \%*\% S \%*\% t(E2)
eps = sig * rnorm(n) # Error
y = matrix.regression.product(X, A) + eps # response
# Estimation
A_sketch <- matrix_sketch(X, y, r, p1, p2)
hatA <- A_sketch[[1]]
norm(hatA - A, type = "F")/norm(A, type = "F") # error
}

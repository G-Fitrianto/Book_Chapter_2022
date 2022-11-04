# Formulation of W Matrix for Lattice Structure (for easier construction)
# Reference: Fitrianto, Tanaka, and Nishii (2018)
# Construct functions for upper and left shifted matrices #

A.mat <- function(dims){
  B <- as(diag(1, dims-1, dims-1), "CsparseMatrix")
  D <- as(matrix(rep(0, dims), nrow = dims, ncol = 1), "CsparseMatrix")
  E <- cbind(rbind(t(D[1:(dims-1),]),B), D)
  return(t(E) + E)
}

library(Matrix); library(MASS)


# Monte Carlo setup
start.time <- Sys.time()
set.seed(11011)
D <- seq(10,150,10)
D.result <- matrix(NA, ncol = 4, nrow = length(D))

for (j in 1:length(D)) {
  r=D[j];c=D[j]
  
  ic <- as(diag(1,ncol = c,nrow = c),"CsparseMatrix")
  ir <- as(diag(1,ncol = r,nrow = r),"CsparseMatrix")
  ac <- A.mat(c); ar <- A.mat(r)
  
  W.s <- kronecker(ar,ic)+kronecker(ir,ac)+kronecker(ar,ac)
  
  rm(ic,ir,ac,ar)
  # Row-standardized W
  
  W <- (1/rowSums(W.s))*W.s
  rm(W.s)
  
  k         <- 100
  delta     <- 0.5
  a         <- 0.8
  n         <- nrow(W)
  In        <- diag(1,nrow = n, ncol = n)
  Inv_W     <- solve(In - delta*W)
  cons      <- matrix(1,ncol = 1,nrow = n)
  
  delta.hat <- numeric(k)
  alpha.hat <- numeric(k)
  
  ## Simulation
  
  for (i in 1:k) {
    e   <- rnorm(n, mean = 0, sd = 1)
    y   <- (a*cons) + Inv_W %*% e
    Wy  <- W %*% y
    
    X   <- cbind(cons, Wy)
    est <- solve(t(X) %*% X) %*% t(X) %*% y
    
    alpha.hat[i] <- est[1,1]
    delta.hat[i] <- est[2,1]
  }
  D.result[j,1] <- a - mean(alpha.hat)
  D.result[j,2] <- delta - mean(delta.hat)
  D.result[j,3] <- var(alpha.hat)
  D.result[j,4] <- var(delta.hat)
}

end.time <- Sys.time()
elapsed <- end.time - start.time
# Bias of OLS

par(mfrow=c(1,1))
plot(D.result[,1], col="blue", lwd=3, main = "Bias of alpha",xlab = "Alpha Parameters",
     ylab = "Magnitude of Bias", type="b", ylim = c(0,max(D.result[,1])))
abline(a=0,b=0,col="black", lty=2, lwd=2)

plot(D.result[,2], col="red", lwd=3, main = "Bias of delta", xlab = "Delta Parameters",
     ylab = "Magnitude of Bias", type="b", ylim = c(min(D.result[,2]),0))
abline(a=0,b=0,col="black", lty=2, lwd=2)

plot(D.result[,3], col="blue", lwd=3, main = "Variance of Alpha Esimation", xlab = "Alpha Parameters",
     ylab = "Variance of alpha.hat", type="b")
abline(a=0,b=0,col="black", lty=2, lwd=2)

plot(D.result[,4], col="red", lwd=3, main = "Variance of Delta Esimation", xlab = "Delta Parameters",
     ylab = "Variance of delta.hat", type="b")
abline(a=0,b=0,col="black", lty=2, lwd=2)

#hist(alpha.hat, freq = F, main="Bias constant parameter")
#hist(delta.hat, freq = F, main="Bias delta parameter")


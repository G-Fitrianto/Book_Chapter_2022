A.mat <- function(dims){
  B <- as.matrix(diag(1, dims-1, dims-1))
  D <- matrix(rep(0, dims), nrow = dims, ncol = 1)
  E <- cbind(rbind(t(D[1:(dims-1),]),B), D)
  return(t(E) + E)
}

library(Matrix)

r=10;c=10

ic <- as.matrix(diag(1,ncol = c,nrow = c))
ir <- as.matrix(diag(1,ncol = r,nrow = r))
ac <- A.mat(c); ar <- A.mat(r)

W.s <- kronecker(ar,ic)+kronecker(ir,ac)+kronecker(ar,ac)

# Row-standardized W

W <- (1/rowSums(W.s))*W.s

k         <- 100
d         <- seq(0,0.999,0.001)
b         <- 0.8
n         <- nrow(W)
In        <- diag(1,nrow = n, ncol = n)

D <- seq(0,0.99,0.01)
D.result <- matrix(NA, ncol = 2, nrow = length(D))

for (j2 in 1:length(D)) {
  delta     <- D[j2]
  Inv_W     <- solve(In - delta*W)
  X         <- matrix(runif(n,0,n*8), ncol=1)
  
  delta.hat  <- numeric(k)
  bias.delta <- numeric(k)
  beta.hat <- numeric(k)
  AIC.temp <- matrix(NA, nrow = k, ncol = length(d))
  BIC.temp <- matrix(NA, nrow = k, ncol = length(d))
  ## Simulation
  
  for (i in 1:k) {
    e   <- rnorm(n, mean = 0, sd = 1)
    y   <- Inv_W %*% (b*X + e)
    Wy  <- W %*% y
    
    for (k1 in 1: length(d)) {
      adjust <- - n*log(det((In - d[k1]*W)))
      By  <- (In - d[k1]*W) %*% y  
      est <- lm(By~0+X)
      
      p    = 1 + 1 + 1
      mse  = sum(est$residuals^2)/ n
      logL = - n * (log(2*pi*exp(1)) + log(mse ))/2
      AIC.temp[i,k1] <- -2*logL + 2*p + adjust
      BIC.temp[i,k1] <- -2*logL + log(n)*p + adjust   
    }
    
    #beta.hat[i]  <- est[1,1]
    d.id <- which(AIC.temp == min(AIC.temp[i,]), arr.ind = T)
    delta.hat[i]  <- d[d.id[2]]
    bias.delta[i] <- delta - d[d.id[2]]
  }
  
  D.result[j2,1] <- delta - mean(delta.hat)
  D.result[j2,2] <- var(delta.hat)
}


#par(mfrow=c(1,2))
#plot(D.result[,1], col="blue", lwd=3, main = "Bias of Delta",xlab = "Delta Parameters",
#     ylab = "Magnitude of Bias", type="b", ylim = c(-4*10^-4,3*10^-4),
#     xaxt="n")
#axis(1, at=seq(1,99,15), labels = seq(0,0.99,0.15))
#abline(a=0,b=0,col="black", lty=2, lwd=2)

#plot(D.result[,2], col="blue", lwd=3, main = "Variance of Delta Estimation", xlab = "Delta Parameters",
#     ylab = "Magnitude of Bias", type="b", ylim = c(0,6*10^-7), xaxt="n")
#axis(1, at=seq(1,99,15), labels = seq(0,0.99,0.15))
#abline(a=0,b=0,col="black", lty=2, lwd=2)

#par(mfrow=c(1,1))
#plot(AIC.temp[12,],lwd=3, xaxt="n", ylab = "Nilai AIC", xlab = "Delta Grid Parameter",
#     main="AIC-OLS Search")
#axis(1, at=seq(1,999,150), labels = seq(0,0.999,0.15))
rm(list = ls())

library(xtable)
library(tinytex)
library(tictoc)
library(doParallel)
registerDoParallel(cores=4)

basedir <- '/home/theodor/Desktop/Courses/ECON_675/Assignments/Assignment\ 2'
setwd(basedir)


######################
##### Question 1 #####
######################

mu1 <- -1.5
sd1 <- sqrt(1.5)
mu2 <- 1
sd2 <- sqrt(1)
mu.true <- 0.5*(mu1+mu2)
sd.true <- sqrt(0.5^2*(sd1^2+sd2^2))

f.true <- function(u) 0.5 * dnorm((u-mu1)/sd1)/sd1 + 0.5 * dnorm((u-mu2)/sd2)/sd2

f.hat <- function(x.eval, X.samp, K, h.n, jack=FALSE){
  if (jack==TRUE) X <- X.samp[X.samp!=x.eval]
  else X <- X.samp
  X.K = K((X-x.eval)/h.n) / h.n
  return(mean(X.K))
}

f.hat.V <- Vectorize(f.hat, vectorize.args="x.eval")

plot.f.hat <- function(X, h.n){
  t <- seq(min(X), max(X), length.out=500)
  y <- numeric(0)
  i <- 0
  for (ti in t){
    i <- i + 1
    y[i] <- f.hat(ti, X, K.epan, h.n)
  }
  plot.new()
  plot(t, y, type="l")
}

# recursive higher order differentiation
DD <- function(expr, name, order = 1) {
  if (order < 1) stop("'order' must be integer >= 1")
  if (order == 1) D(expr, name)
  else DD(D(expr, name), name, order - 1)
}

# evaluation of symbolic derivative
DDeval <- function(fxn, u, o=1) eval(DD(fxn, "u", order=o))

# Epanechnikov kernel
K.epan <- function(u) (0.75 * (1-u^2) * (abs(u)<=1))

mu.l <- function(l, fxn){
  mu.l.int <- function(u) (u^l * fxn(u))
  return(integrate(mu.l.int, -Inf, Inf))
}

theta.l <- function(l, fxn){
  if (l==0) {
    theta.l.int <- function(u) (fxn(u))^2
  }
  else {
    # CAREFUL: symbolic derivatives do not work on K because of indicator function
    # should add option for numerical derivative
    theta.l.int <- function(u) (DDeval(body(fxn), u, o=l))^2
  }
  return(integrate(theta.l.int, -Inf, Inf))
}

# AIMSE optimal bandwidth
h.AIMSE <- function(s, P, K, f, n){
  return( ( (2*s + 1) * (factorial(P)^2) / (2*P) * theta.l(s, K)[[1]] / ( (theta.l(s+P, f)[[1]]) * (mu.l(P, K)[[1]])^2 ) * 1 / n )^(1 / (2*s + 2*P + 1)) )
}

h.hat.AIMSE <- function(s, P, K, f, n){
  return( ( theta.l(s, K)[[1]] / ( (theta.l(s+P, f)[[1]]) * (mu.l(P, K)[[1]])^2 ) * 1 / n )^(1 / (2*s + 2*P + 1)) )
}

draw <- function(n) {
  X <- rep(NaN, n)
  N1 <- rbinom(1, size=n, prob=0.5)
  X[1:N1] <- rnorm(N1, mean=mu1, sd=sd1)
  X[(N1+1):n] <- rnorm(n-N1, mean=mu2, sd=sd2)
  return(X)
}

##### Question 1

#set.seed(90210)

n <- 1000
s <- 0
P <- 2

### 1.3.a

h.opt <- h.AIMSE(s, P, K.epan, f.true, n)

### 1.3.b

M <- 1000
H <- h.opt * seq(0.5, 1.5, 0.1)
nh <- length(H)

tic("Data simulation")
MSE.L1 <- matrix(NaN, M, nh)
MSE.L0 <- matrix(NaN, M, nh)
h.hat <- numeric(0)
cat('Simulating data...\n')
for(m in 1:M) {
  if (m%%100 == 0) cat('---Draw ', format(m), ' out of ', format(M), '\n')
  X <- draw(n)
  
  ## calculate MSE
  
  # run sequentially for each bandwidth
  # j <- 0
  # for(h in H) {
  #   j <- j + 1
  #   MSE.L1[m, j] <- mean((f.hat.V(X, X, K.epan, h) - f.true(X))^2)
  #   MSE.L0[m, j] <- mean((f.hat.V(X, X, K.epan, h, jack=TRUE) - f.true(X))^2)
  # }
  
  # run in parallel
  results = foreach(j=1:nh, .combine=data.frame) %dopar% {
    c(mean((f.hat.V(X, X, K.epan, H[j]) - f.true(X))^2), mean((f.hat.V(X, X, K.epan, H[j], jack=TRUE) - f.true(X))^2))
  }
  MSE.L1[m, ] <- as.numeric(results[1, ])
  MSE.L0[m, ] <- as.numeric(results[2, ])
  
  ### 1.3.d
  # calculate
  mu.est <- mean(X)
  sd.est <- sd(X)
  f.est <- function(u) dnorm((u-mu.est)/sd.est)/sd.est
  h.hat[m] <- h.hat.AIMSE(s, P, K.epan, f.est, n)
}
IMSE.L1 <- colMeans(MSE.L1)
IMSE.L0 <- colMeans(MSE.L0)
h.min.L1 <- H[which.min(IMSE.L1)]
h.min.L0 <- H[which.min(IMSE.L0)]
h.bar <- mean(h.hat)
toc()

# create figure
plot(H, IMSE.L1, type="l", xlab="bandwidth", ylab="IMSE")
lines(H, IMSE.L0, type="l", lty=2)
lines(rep(h.opt, 2), c(-1, 1), type="l", col="red", lwd=3)
lines(rep(h.min.L1, 2), c(-1, 1), type="l", lty=1)
lines(rep(h.min.L0, 2), c(-1, 1), type="l", lty=2)
lines(rep(h.bar, 2), c(-1, 1), type="l", col="orange", lwd=2)


######################
##### Question 2 #####
######################

# generate the polynomial basis
gen.P = function(Z,K) {
  if (K==0)   out = NULL;
  if (K==1)   out = poly(Z,degree=1,raw=TRUE);
  if (K==2)  {out = poly(Z,degree=1,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^2);}
  if (K==2.5) out = poly(Z,degree=2,raw=TRUE);
  if (K==3)  {out = poly(Z,degree=2,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^3);}
  if (K==3.5) out = poly(Z,degree=3,raw=TRUE);
  if (K==4)  {out = poly(Z,degree=3,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^4);}
  if (K==4.5) out = poly(Z,degree=4,raw=TRUE);
  if (K==5)  {out = poly(Z,degree=4,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^5);}
  if (K==5.5) out = poly(Z,degree=5,raw=TRUE);
  if (K>=6)  {out = poly(Z,degree=5,raw=TRUE); for (k in 6:K) for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^k);}
  ## RETURN POLYNOMIAL BASIS
  return(out)
}

# DGP
data_gen <- function(n) {
  X <- matrix((runif(n*5)-0.5)*2, ncol=5)
  V <- matrix(rnorm(n), ncol=1)
  G <- matrix(exp(diag(X %*% t(X))), ncol=1)
  E <- matrix(0.3637899 * (1 + diag(X %*% t(X))) * V, ncol=1)
  U <- matrix(rnorm(n), ncol=1)
  tt <- matrix(sqrt(diag(X %*% t(X))) + U >= 0, ncol=1) * 1
  Y <- matrix(tt + G + E, ncol=1)
  return(list(Y=Y, X=X, tt=tt))
} 

### 2.5.a

n   <- 500
K   <- c(1, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10)
K.r <- c(6, 11, 21, 26, 56, 61, 126, 131, 252, 257, 262, 267, 272, 277)
nK  <- length(K)
M   <- 1000
theta.hat <- matrix(NaN, ncol=nK, nrow=M)
se.hat    <- theta.hat

### 2.5.a
set.seed(90210)
tic("Simulation 2.5.b")
for (m in 1:M) {
  if (m%%100 == 0) cat('---Draw ', format(m), ' out of ', format(M), '\n')
  data <- data_gen(n)
  X <- data$X
  Y <- data$Y
  tt <- data$tt
  for (k in 1:nK) {
    X.pol <- cbind(1, gen.P(X, K[k]))
    X.Q   <- qr.Q(qr(X.pol))
    MP     <- diag(rep(1,n)) - X.Q %*% t(X.Q)
    Y.M <- MP %*% Y
    tt.M <- MP %*% tt
    theta.hat[m, k] <- (t(tt.M) %*% Y.M) / (t(tt.M) %*% tt.M)
    Sigma <- diag((as.numeric((Y.M - tt.M*theta.hat[m, k])))^2)
    se.hat[m, k] <- sqrt(t(tt.M) %*% Sigma %*% tt.M) / (t(tt.M) %*% tt.M)
  }
}
toc()

table <- matrix(NaN, ncol=6, nrow=nK)
for (k in 1:nK) {
  table[k, 1] <- K.r[k]
  table[k, 2] <- mean(theta.hat[, k]) - 1                           # bias
  table[k, 3] <- sd(theta.hat[, k])                                 # standard deviation
  table[k, 4] <- table[k, 2]^2 + table[k, 3]^2                      # mse
  table[k, 5] <- mean(se.hat[, k])                                  # mean standard error
  table[k, 6] <- mean((theta.hat[, k] - 1.96 * se.hat[, k] > 1) | 
                        (theta.hat[, k] + 1.96 * se.hat[, k] < 1))  # rejection rate
}
write.table(round(table,3), "partial_linear.txt", sep=" & ", eol="\\\\ \n", col.names = FALSE, row.names = FALSE)

### 2.5.c

n   <- 500
K   <- c(1, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10)
K.r <- c(6, 11, 21, 26, 56, 61, 126, 131, 252, 257, 262, 267, 272, 277)
nK  <- length(K)
M   <- 1000

# generate the polynomial basis
gen.P = function(Z,K) {
  if (K==0)   out = NULL;
  if (K==1)   out = poly(Z,degree=1,raw=TRUE);
  if (K==2)  {out = poly(Z,degree=1,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^2);}
  if (K==2.5) out = poly(Z,degree=2,raw=TRUE);
  if (K==3)  {out = poly(Z,degree=2,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^3);}
  if (K==3.5) out = poly(Z,degree=3,raw=TRUE);
  if (K==4)  {out = poly(Z,degree=3,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^4);}
  if (K==4.5) out = poly(Z,degree=4,raw=TRUE);
  if (K==5)  {out = poly(Z,degree=4,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^5);}
  if (K==5.5) out = poly(Z,degree=5,raw=TRUE);
  if (K>=6)  {out = poly(Z,degree=5,raw=TRUE); for (k in 6:K) for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^k);}
  ## RETURN POLYNOMIAL BASIS
  return(out)
}

# DGP
data_gen <- function(n) {
  X <- matrix((runif(n*5)-0.5)*2, ncol=5)
  V <- matrix(rnorm(n), ncol=1)
  G <- matrix(exp(diag(X %*% t(X))), ncol=1)
  E <- matrix(0.3637899 * (1 + diag(X %*% t(X))) * V, ncol=1)
  U <- matrix(rnorm(n), ncol=1)
  tt <- matrix(sqrt(diag(X %*% t(X))) + U >= 0, ncol=1) * 1
  Y <- matrix(tt + G + E, ncol=1)
  return(list(Y=Y, X=X, tt=tt))
} 

# cross validation function
K.CV <- function(tt, X, Y) {
  temp <- rep(NaN, nK)
  for (k in 1:nK) {
    X.pol <- cbind(1, tt, gen.P(X, K[k]))
    X.Q   <- qr.Q(qr(X.pol))
    XX <- X.Q %*% t(X.Q)
    Y.hat <- XX %*% Y
    W <- diag(XX)
    temp[k] <- mean(((Y-Y.hat) / (1-W))^2)
  }
  return(which.min(temp))
}

theta.hat2 <- rep(NaN, M)
se.hat2    <- theta.hat2
K.hat2     <- theta.hat2

set.seed(90210)
tic("Simulation 2.5.c")
for (m in 1:M) {
  if (m%%100 == 0) cat('---Draw ', format(m), ' out of ', format(M), '\n')
  data <- data_gen(n)
  X <- data$X; Y <- data$Y; tt <- data$tt
  k.opt <- K.CV(tt, X, Y)
  X.pol <- cbind(1, gen.P(X, K[k.opt]))
  X.Q   <- qr.Q(qr(X.pol))
  MP    <- diag(rep(1,n)) - X.Q %*% t(X.Q)
  Y.M   <- MP %*% Y
  tt.M  <- MP %*% tt
  theta.hat2[m] <- (t(tt.M) %*% Y.M) / (t(tt.M) %*% tt.M)
  Sigma         <- diag((as.numeric((Y.M - tt.M*theta.hat[m, k])))^2)
  se.hat2[m]    <- sqrt(t(tt.M) %*% Sigma %*% tt.M) / (t(tt.M) %*% tt.M)
  K.hat2[m]     <- K.r[k.opt]
}
proc.time() - ptm

# summary of the cross validation 
table(K.hat2)
# estimator
summary(theta.hat2)
sd(theta.hat2)
summary(se.hat2)
sd(se.hat2)

par(mfrow=c(1,2))
hist(theta.hat2, freq=FALSE, xlab="theta-hat", ylab="", main="")
lines(c(mean(theta.hat2), mean(theta.hat2)), c(-1, 20), col="red", lwd=3)
hist(se.hat2, freq=FALSE, xlab="s.e.", ylab="", main="")
lines(c(mean(se.hat2), mean(se.hat2)), c(-1, 80), col="red", lwd=3)

par(mfrow=c(1,2))
CI.l <- theta.hat2 - 1.96 * se.hat2
CI.r <- theta.hat2 + 1.96 * se.hat2
# rejection rate
mean(1 < CI.l | 1 > CI.r)
plot(1:M, CI.l, type="l", ylim=c(0,2), xlab="simulations", ylab="CI")
lines(1:M, CI.r)
abline(1, 0, col="red", lwd=2)

temp <- sort(CI.l, index.return=TRUE)
CI.l <- temp$x
CI.r <- CI.r[temp$ix]
plot(1:M, CI.l, type="l", ylim=c(0,2), xlab="simulations", ylab="CI")
lines(1:M, CI.r)
abline(1, 0, col="red", lwd=2)

### Repeat 2.5.c for derivative of mu


######################
##### Question 3 #####
######################

#rm(list=ls())

##### 3.4

# Setup
theta <- 1    # true value of parameter of interest
d     <- 5    # x_i in R^d
n     <- 500  # observations
M     <- 1000 # replications

# True function
L2norm <- function(x) rowSums(x^2)
g0 <- function(x) exp(L2norm(x))

# Constructing basis
r <- function(K){
  if(K==6){return(r6)}
  if(K==11){return(r11)}   
  if(K==21){return(r21)}   
  if(K==26){return(r26)}  
  if(K==56){return(r56)}   
  if(K==61){return(r61)}   
  if(K==126){return(r126)}   
  if(K==131){return(r131)}   
  if(K==252){return(r252)}   
  if(K==257){return(r257)}   
  if(K==262){return(r262)}   
  if(K==267){return(r267)}   
  if(K==272){return(r272)}   
  if(K==277){return(r277)}   
}

# Calculations
est <- function(y, t, K){
  e     <- rep(0,K+1)
  e[1]  <- 1 
  rmat  <- r(K) 
  Z     <- cbind(t, rmat)
  
  # (i) thetahat (OLS)
  M1 <- t(Z) %*% Z
  M2 <- t(Z) %*% y
  M3 <- solve(M1)
  betahat  <- M3 %*% M2
  thetahat <- e %*% betahat
  
  # (ii) bias
  bias  <- thetahat - 1
  
  # (iv) Vhat
  Mp    <- diag(length(y)) - rmat%*%solve(t(rmat)%*%rmat)%*%t(rmat)
  epshat <- (y - Z%*%betahat)
  sigmahat <- diag(array(epshat^2))
  Vhat <- solve(t(t)%*%Mp%*%t)%*%(t(t)%*%Mp%*%sigmahat%*%Mp%*%t)%*%solve(t(t)%*%Mp%*%t)
  
  # (v) coverage
  CIl   <- thetahat + sqrt(Vhat)*qnorm(.025)
  CIu   <- thetahat + sqrt(Vhat)*qnorm(.975)
  cover <- ((CIl <= 1) & (CIu >= 1))*1
  
  # extra: MSE
  W     <- diag(Z%*%(M3)%*%t(Z))
  yhat  <- Z%*%betahat
  mse   <- mean(((y-yhat)/(1-W))^2)
  
  # return a row vector with thetahat, bias, vhat, cover for given K
  return(cbind(thetahat, bias, Vhat, cover, mse))
}

### 3.4.a

# Polynomial basis expansions to calculate K = ...
kvals <- c(6, 11, 21, 26, 56, 61, 126, 131, 252, 257, 262, 267, 272, 277)

# Matrices to fill for K estimations over M iterations
thetahat  <- matrix(NA, ncol = length(kvals), nrow = M)
bias      <- matrix(NA, ncol = length(kvals), nrow = M)
Vhat      <- matrix(NA, ncol = length(kvals), nrow = M)
cover     <- matrix(NA, ncol = length(kvals), nrow = M)
mse       <- matrix(NA, ncol = length(kvals), nrow = M)
kstar     <- rep(NA,M)

# Simulation
tic("Simulation 3.4")
for(m in 1:M){
  if (m%%100 == 0) cat('---Draw ', format(m), ' out of ', format(M), '\n')
  # Generate data
  x   <- matrix(runif(n=n*d,min = -1, max = 1), nrow = n, ncol = d)
  v   <- rnorm(n)
  u   <- rnorm(n)
  eps <- 0.3637899*(1 + L2norm(x))*v
  g   <- g0(x)
  t   <- matrix((sqrt(L2norm(x)) + u >=0)*1, ncol=1)
  y   <- t*theta + g + eps
  
  r6   <- cbind(rep(1,n),x)
  r11  <- cbind(r6,x^2)
  r21  <- cbind(rep(1,n),poly(x,degree = 2, raw = T))
  r26  <- cbind(r21,x^3)
  r56  <- cbind(rep(1,n),poly(x,degree = 3, raw = T))
  r61  <- cbind(r56,x^4)
  r126 <- cbind(rep(1,n),poly(x,degree = 4, raw = T))
  r131 <- cbind(r126,x^5)
  r252 <- cbind(rep(1,n),poly(x,degree = 5, raw = T))
  r257 <- cbind(r252,x^6)
  r262 <- cbind(r257,x^7)
  r267 <- cbind(r262,x^8)
  r272 <- cbind(r267,x^9)
  r277 <- cbind(r272,x^10)
  
  ### 3.4.b
  
  # (i)   avg of thetahat(K)
  # (ii)  avg bias of thetahat(K)
  # (iii) sample variance of thetahat(K)
  # (iv)  avg of Vhat_HCO 
  # (v)   average coverate rate of the confidence intervals 
  
  temp <- matrix(NA,nrow = length(kvals),ncol = 4)
  for(k in 1:length(kvals)){
    knum <- as.numeric(kvals[k])
    temp <- est(y, t, knum)
    thetahat[m, k] <- temp[1]
    bias[m, k] <- temp[2]
    Vhat[m, k] <- temp[3]
    cover[m, k]<- temp[4]
    mse[m, k]  <- temp[5]
  }
  
  kstar[m] <- kvals[which.min(mse[m,])]
}
toc()

avg.thetahat  <- colMeans(thetahat)  
avg.bias      <- colMeans(bias)
avg.var       <- apply(thetahat, 2, sd)
avg.Vhat      <- colMeans(Vhat)
avg.cover     <- colMeans(cover)
avg.kstar     <- mean(kstar)

### 3.4.c

thetahat_CV   <- rep(NA,M)
bias_CV       <- rep(NA,M)
var_CV        <- rep(NA,M)
Vhat_CV       <- rep(NA,M)
cover_CV      <- rep(NA,M)

for(m in 1:M){
  thetahat_CV[m] <- thetahat[m,which.min(mse[m,])]
  bias_CV[m] <- bias[m,which.min(mse[m,])]
  Vhat_CV[m] <- Vhat[m,which.min(mse[m,])]
  cover_CV[m] <- cover[m,which.min(mse[m,])]
}

thetahat_KCV  <- mean(thetahat_CV)
bias_KCV      <- mean(bias_CV)
var_KCV       <- sd(thetahat_CV)
Vhat_KCV      <- mean(Vhat_CV)
cover_KCV     <- mean(cover_CV)

### Output

table <- cbind(kvals,avg.thetahat,avg.bias,avg.var,avg.Vhat,avg.cover)
colnames(table) <- c("K","thetahat","bias","var","Vhat","coverage")
table3.4.b <- xtable(table, caption = "Results by K value")
print(table3.4.b,type = "latex", file = "solutions/table_3.4.b.tex")

table <- cbind(avg.kstar, thetahat_KCV, bias_KCV, var_KCV, Vhat_KCV, cover_KCV)
colnames(table) <- c("K_CV","thetahat","bias","var","Vhat","coverage")
table3.4.c <- xtable(table, caption = "Results using CV approach")
print(table3.4.c,type = "latex", file = "solutions/table_3.4.c.tex")

table <- rbind(cbind(avg.thetahat,avg.bias,avg.var,avg.Vhat,avg.cover),cbind(thetahat_KCV, bias_KCV, var_KCV, Vhat_KCV, cover_KCV))
rownames(table) <- c(kvals,round(avg.kstar,1))
colnames(table) <- c("thetahat","bias","var","Vhat","coverage")
table3.4 <- xtable(table, caption = "Results by K value with $K_CV$ last")
print(table3.4,type = "latex", file = "solutions/table_3.4cR.tex")





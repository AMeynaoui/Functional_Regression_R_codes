##                    Functions for the simulation of samples (i), (ii) and (iii)                            ##


simuX=function(n,p=100,kreal=50,lambda=(1:kreal)^(-2),set="CM13",loi="norm",...)
{
  
  # generate a n*p matrix (X_i(t_j)) i=1,...,n; j=1,...,p from a truncated Karhunen-Loeve decomposition
  # X_i(t_j)=sum_{k=1}^{k_real} ksi[k] sqrt(lambda[k]) phi_k(t_j)
  # where 
  ## {ksi_k}_{k=1}~N(0,1) i.i.d., if noise = "norm"
  ## {ksi_k}_{k=1}~U(-sqrt(3),sqrt(3))if noise = "unif"
  # lambda is a vector of size k_real positive numbers 
  # phi_k are functions and t_1,...,t_p is a regular grid of [0,1]
  ## phi_k(t) = sqrt(2)*sin((j-0.5)pi t) if set="CM13"
  ## phi_k(t) = sqrt(2)*cos(j pi t) if set = "IK18"
  
  if (loi=='norm')  ksi=rnorm(n*kreal) else ksi=runif(n*kreal,-sqrt(3),sqrt(3))
  ksi=matrix(ksi,n,kreal)*matrix(rep(sqrt(lambda),n),n,kreal,byrow=TRUE)
  x=matrix(rep(seq(0,1,length.out=p),kreal),kreal,p,byrow=TRUE)
  if (set == "CM13"){
    Jm=matrix(rep(seq(0.5,kreal-0.5,by=1),p),kreal,p) 
    base=sqrt(2)*sin(pi*Jm*x)
  } else {
    Jm=matrix(rep(seq(1,kreal,by=1),p),kreal,p) 
    base=sqrt(2)*cos(pi*Jm*x)
  }
  ksi%*%base 
}

MvB=function(n,p)
{
  X=rnorm(n*p)/sqrt(p)
  X=matrix(X,ncol=p)
  MvB=matrix(0,ncol=p,nrow=n)
  for (i in 1:n)
  {
    MvB[i,]=cumsum(X[i,])
  }
  MvB
}

simuY <- function(X,S,noise="Brownian",q,fac=1,delta=1.1,k1=50,...){
  # generate a n*q matrix (X_i(t_j)) i=1,...,n; j=1,...,q such that
  # Y_i(t_j) \approx S X_i (t_j)
  # where
  ## X is a n*p matrix containing the X_i's
  ## S is an operator (i.e. a function of a function f and t)
  ## noise is a brownian motion B_t/fac if noise = "Brownian"
  ## noise = sum_{k=1}^{k_1} k^{-delta/2}\xi_k\phi_k(t_j) with \phi_k(t_j)=sqrt(2) cos(2 pi k t_j) if noise = "IK18"
  
  discrY <- seq(0,1,length.out = q) # grid for the Yi's
  MatSX <- S(discrY,X,...) #  matrix of the SXi's
  if (noise=="Brownian"){
    epsilon = MvB(n,q)/fac
  } else {
    epsilon = simuX(n,q,lambda = (1:k1)^(-delta/2),set="CM13",noise="unif")
  }
 MatSX + epsilon
}

operatorCM13 <- function(t,X){ # operator S (t vector of points of a regular grid) and X a n*p matrix of n discretized function)
  p = ncol(X)
  q = length(t)
  mats = matrix(rep(seq(0,1,length.out=p),q),p,q)
  matt = matrix(rep(t,p),p,q,byrow = T )
  X %*% kernel(mats,matt)/p
}

operatorIK18 <- function(t,X,beta=3,gam=2.5){ # operator S (t vector of points of a regular grid) and X a n*p matrix of n discretized function)
  p = ncol(X)
  q = length(t)
  gridX = seq(0,1,length.out=p)
  Phip = sqrt(2)*cos(2*pi*tcrossprod(1:k1,gridX))
  B = 4*tcrossprod((-1)^(1:k1)*(1:k1)^(-gam),(-1)^(1:k1)*(1:k1)^(-beta))
  Phiq = sqrt(2)*cos(2*pi*tcrossprod(1:k1,t))
  crossprod(tcrossprod(Phip,X),crossprod(B,Phiq))
}
##                    Functions for the estimation of the slope operator S                                  ##

norm2 <- function(x){sqrt(sum(x^2))} # usual norm of R^length(x)

Estim_Sn <- function(MatX,MatY,k1){ # Estimation of S
  
  # Entries : 
  ## MatX : n*p matrix containing (X_i(t_j)), i=1,...,n; j=1,...,p (t_j)_j a regular grid
  ## MatY : n*q matrix containing (Y_i(s_j)), i=1,...,n; j=1,...,q (s_j)_j a regular grid
  ## k1 : dimension of the projection space
  
  # Output :
  ## matrix of the estimated kernel of S (Sronde)
  
  p <- ncol(MatX) # number of discretization points of X
  q <- ncol(MatY) # number of discretization points of Y
  n <- nrow(MatX) # number of individuals
  MatSn <- matrix(NA, nrow = q, ncol = p) 
  
  PCAX <- fPCA(MatX)  # fPCA basis of the X_i's
  Eigen_val <-  PCAX$hlambda # vector of eigenvalues
  Eigen_fun <-  PCAX$hPsi    # matrix of discretized eigenfunctions 
  
  Xe <- PCAX$coeffs[,1:k1]   # score matrix (coefficients of the X_i's in the fPCA basis)
  Xe_lambda <- sweep(Xe,2,Eigen_val[1:k1],"/")  # matrix of the coefficients divided by the estimated eigenfunctions
  
  for (i in 1:q){
    for(j in 1:p){
      A <- sweep(Xe_lambda,1,MatY[,i],"*")
      A <- sweep(A,2,Eigen_fun[1:k1,j],"*")
      MatSn[i,j] <- 1/n*sum(A)
    }
  }
  return(MatSn)
}

Dim.Selec <- function(kappa,Nn,MatX,MatY,sigma2){# dimension selection function whith known noise variance
  # Entries :
  ##kappa : calibration parameter 
  ##Nn : maximal dimension
  ## MatX : n*p matrix containing (X_i(t_j)), i=1,...,n; j=1,...,p (t_j)_j a regular grid
  ## MatY : n*q matrix containing (Y_i(s_j)), i=1,...,n; j=1,...,q (s_j)_j a regular grid
  ## sigma2 : noise variance 
  
  # Output : a list a two 
  ## $Dim_Optim : selected dimension
  ## $Risk.Em : value of the criterion (contrast + penalty) for the selected dimension (minimum of the criterion)
  
  p <- ncol(MatX) # number of discretization points of X
  q <- ncol(MatY) # number of discretization points of Y
  n <- nrow(MatX) # number of individuals
  Risk.Em <- rep(NA,Nn-1) # vector of contrast values (as a function of dimension)
  for(m in 2:Nn){
    MatSm <- Estim_Sn(MatX,MatY,m) # estimator on the m-first principal components
    Shift.Y.SX <- MatY - 1/p*t(MatSm%*%t(MatX)) # matrix of Yi-SXi (residuals)
    Mat.Gn.m <- apply(Shift.Y.SX,1,norm2) # matrix of the squared norms of (Yi - SXi)
    Gn_Sm <- 1/n*sum(Mat.Gn.m) # value of the contrast Gamma_n(Sm)
    Risk.Em[m-1] <- Gn_Sm + kappa*sigma2*m/n # value of the criterion for the dimension m
  }
  Dim_Optim <- as.numeric(which(Risk.Em == min(Risk.Em),arr.ind = T)) + 1  # optimal dimension
  return(list(Dim.Optim=Dim_Optim,Risk.Em=Risk.Em))
}


Dim.Selec.sigmahat <- function(kappa,Nn,MatX,MatY){ # dimension selection function when the noise variance is not known
  # Entries :
  ##kappa : calibration parameter 
  ##Nn : maximal dimension
  ## MatX : n*p matrix containing (X_i(t_j)), i=1,...,n; j=1,...,p (t_j)_j a regular grid
  ## MatY : n*q matrix containing (Y_i(s_j)), i=1,...,n; j=1,...,q (s_j)_j a regular grid
  
  # Output : a list a two 
  ## $Dim_Optim : selected dimension
  ## $Risk.Em : value of the criterion (contrast + penalty) for the selected dimension (minimum of the criterion)
  
  
  p <- ncol(MatX) # number of discretization points of X
  q <- ncol(MatY) # number of discretization points of Y
  n <- nrow(MatX) # number of individuals
  
  Risk.Em <- rep(NA,Nn-1) # vector of contrast values (as a function of dimension)
  for(m in 2:Nn){
    MatSm <- Estim_Sn(MatX,MatY,m) # estimator on the m-first principal components
    Shift.Y.SX <- MatY - 1/p*t(MatSm%*%t(MatX)) # matrix of Yi-SXi (residuals)
    Mat.Gn.m <- apply(Shift.Y.SX,1,norm2) # matrix of the squared norms of (Yi - SXi)
    Gn_Sm <- 1/n*sum(Mat.Gn.m) # value of the contrast Gamma_n(Sm)
    Risk.Em[m-1] <- Gn_Sm + kappa*Gn_Sm*m/n # value of the criterion for the dimension m
  }
  Dim_Optim <- as.numeric(which(Risk.Em == min(Risk.Em),arr.ind = T)) + 1  # optimal dimension
  return(list(Dim.Optim=Dim_Optim,Risk.Em=Risk.Em))
}

CVpred <- function(X,Y,kappa,Nn){
  # Entries :
  ## X: n*p matrix containing (X_i(t_j)), i=1,...,n; j=1,...,p (t_j)_j a regular grid
  ## Y: n*q matrix containing (Y_i(s_j)), i=1,...,n; j=1,...,q (s_j)_j a regular grid
  ## kappa: calibration parameter 
  ## Nn: maximal dimension
  
  
  # Output: a list of three
  ## $pred_error: the cross-validated prediction error
  ## $dim_selec: dimension selected with the sample {(X_j,Y_j)} j!=i
  ## $Yhat: n*q matrix with the cross-validated predictions of the Yi's
  
  
  n = nrow(X)
  p = ncol(X)
  
  mhat = rep(NA,n)
  risk = rep(NA,n)
  Yhat = matrix(NA,n,p)
  ntest = 1
  
  for (i in 1:n){
    napp = n-ntest
    Xapp = X[-i,]
    Yapp = Y[-i,]
    Xtest = X[i,]
    Ytest = Y[i,]
    Selecapp <- Dim.Selec.sigmahat(kappa,Nn,Xapp,Yapp)
    mhat[i]=Selecapp$Dim.Optim
    hatSn = Estim_Sn(Xapp,Yapp,Selecapp$Dim.Optim)
    Yhat[i,]=1/p*hatSn%*%Xtest
  }
  
  list(pred_error = sqrt(rowMeans((Yhat-Y)^2)),dim_selec=mhat,Yhat=Yhat)
}
#                           Simulation study of the article                               #
# Adaptive nonparametric estimation in the functional linear model with functional output #
# by Gaelle Chagny, Anouar Meynaoui and Angelina Roche                                    #

# Libraries and functions

library("rstudioapi")
setwd(dirname(rstudioapi::callFun("getActiveDocumentContext")$path))
source("PCA.R") # loading the functions for fPCA
source("functions_simusXY.R")
source("functions_estim.R")

# Parameters

n <- 1000 # size of sample
p <- 100 # size of grid for the Xi's
q <- 100 # size of grid for the Yi's

# Data generation 

## Setting (i)

Xi = simuX(n,kreal=8,lambda = (pi*(0.5:7.5))^(-2))
kernel <- function(s,t) {s^2 + t^2} # kernel Srond of operator S

Yi = simuY(Xi,operatorCM13,noise="Brownian",q=q,fac=20)

## Setting (ii)

Xii = Xi
Yii = simuY(Xii,operatorCM13,noise="Brownian",q=q,fac=2)

## Setting (iii)
alpha = 1.2
beta = 3
gam = 2.5
delta = 1.1
k1 = 50

Xiii = simuX(n,lambda = (1:k1)^(-alpha/2),set="IK18",loi="unif")
Yiii = simuY(Xiii,operatorIK18,noise="IK18",q=q,fac=2,beta=beta,gam=gam,delta=delta)

# Estimation (including dimension selection)

kappa = 0.6
Nn = 40
sigma2i = 1/(2*(20)^2)
sigma2ii = 1/8
sigma2iii = sum((1:k1)^(-delta))
hatm1i <- Dim.Selec(kappa,50,Xi,Yi,sigma2i)
Sni <- Estim_Sn(Xi,Yi,hatm1i$Dim.Optim)


fPCA <- function(X){
# prend en entree une matrice n lignes, p colonnes des donnees discretisees
# retourne une liste contenant : 
  # - une matrice a n-1 lignes, p colonnes des estimateurs de la base des fonctions propres
  # - un vecteur contenant les valeurs propres. 
  # - une matrice n lignes, p colonnes des coefficients des X_i dans la base de l'ACP
  Xc = t(t(X) - rowMeans(t(X)))
  p<-dim(X)[2]
  n<-dim(X)[1]
  coeffs<-Xc/sqrt(p) # coefficients de X_i-E[X_i] dans la base des histogrammes à p bins phi_j(t)=sqrt(p)1_[(j-1)/p,j/p]
  Gn<-crossprod(coeffs)/n # approximation de la matrice <\hat\Gamma phi_j,phi_k>
  ev<-eigen(Gn,symmetric=TRUE)
  vecrenorm = sqrt(p)*t(ev$vec) # les vecteurs propres de Gn nous donnent les vecteurs <\hat e_j,\phi_k>
  # donc \hat e_j(t_k)\approx \sqrt(p)<\hat e_j,\phi_k>
  # coefficients de X dans la base hPsi
  coeffs = coeffs%*%ev$vec
  list(hlambda=ev$val,hPsi=vecrenorm,coeffs=coeffs)
}


MvB=function(n,p)
{
  # Retourne une matrice de taille n*p (X_i(t_j)) i=1,...,n, j=1,...,p
  # tel que X_i est un mouvement brownien standard.
  X=rnorm(n*p)/sqrt(p)
  X=matrix(X,ncol=p)
  MvB=matrix(0,ncol=p,nrow=n)
  for (i in 1:n)
  {
    MvB[i,]=cumsum(X[i,])
  }
  MvB
}

# n=500
# p=100
# X = MvB(n,p)
# t = seq(0,1,length.out=p)
# 
# resPCA = fPCA(X)
# # Valeurs propres 
# resPCA$hlambda
# (pi*seq(0.5,100.5,by=1))^(-2)
# 
# # Fonctions propres
# par(mfrow=c(3,3))
# for (j in 1:9){
#   plot(t,resPCA$hPsi[j,],type='l',main=paste(j,"eme fonction propre"),xlab='t',ylab='Psi_j(t)')
#   points(t,sign(resPCA$hPsi[j,1])*sqrt(2)*sin(pi*(j-0.5)*t),type='l',col=2)
# }
# 
# # Coefficients 
# par(mfrow=c(1,1))
# plot(resPCA$coeffs[,1],resPCA$coeffs[,2],xlab='<X_i,psi_1>',ylab='<X_i,psi_2>',main='Representation des individus sur les deux premieres CP')
# 
# var(resPCA$coeffs[,1])
# resPCA$hlambda[1]

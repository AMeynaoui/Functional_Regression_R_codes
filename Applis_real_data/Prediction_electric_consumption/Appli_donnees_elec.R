rm(list=ls())
source("../../PCA.R") 
source("../../functions_estim.R")


# Data
energy_data <- read.csv("energydata_complete.csv")
View(energy_data)

energy_data$date <- strptime(as.character(energy_data$date),format="%Y-%m-%d %H:%M:%S")
energy_data$date <- as.POSIXct(energy_data$date,tz = "UTC")

dim(energy_data)
names(energy_data)

n = 137
t = seq(from=0,by=24/144,length=144)

# From the temporal series energy_data[,2] to a matrix of functional data
energy_data_day = matrix(energy_data[1:19728,2],137,144,byrow=TRUE)


# Covariable X: log consumption of day j recentered
Xinit = energy_data_day[1:(n-1),]
Xinitbar = colMeans(Xinit)
Xbar = log(colMeans(Xinit))
X = log(Xinit)-Xbar

plot(range(t),range(Xinit),type='n',xlab='time(hour)',ylab='Electric consumption (kWh)')
for (i in 1:(n-1)){
  points(t,Xinit[i,],type='l',col='gray')
}
points(t,Xinitbar,type='l',col='red',lty=2)

plot(range(t),range(X),type='n',xlab='time(hour)',ylab=expression(X[i]))
for (i in 1:(n-1)){
  points(t,X[i,],type='l',col='lightblue')
}
abline(h=0,col='red',lty=2)

# Variable Y: log consumption of day j+1 recentered 
Yinit = energy_data_day[2:n,]
Ybar = colMeans(log(Yinit))
Y = log(Yinit)-Ybar
n = nrow(X)
p = ncol(X)


kappa = 0.6
Nn = 15
# Cross-validation
kappa = 0.6
Nn = 15


res.CV <- CVpred(X,Y,kappa,Nn)
par(mfrow=c(1,1))
barplot(prop.table(table(res.CV$dim_selec)),col='blue',xlab="Dimension selected",ylab='Frequency')
hist(res.CV$pred_error,freq=F,col='orange',xlab='L^2 prediction error',main='')


par(mfrow=c(1,3))
ord = sort.int(pred_error,index.return = T)$ix
show = ord[c(1,round(n/2),n)]
par(mfrow=c(1,3))
for (i in show){
  plot(t,Y[i,],type='l',xlab='t',ylab='Prediction of Y',col='lightblue',main=paste("i=",i),lwd=2)
  points(t,res.CV$Yhat[i,],type='l',col='darkgreen',lty=2,lwd=2)
}
legend('topleft',c(expression(Y[i]),expression(hat(Y)[i]^(-i))),lty=c(1,2),lwd=c(2,2),col=c('lightblue','darkgreen'))

# retour a la consommation (decentrage et passage a l'exponentielle)

for (i in show){
  Yhattestrec = exp(Yhat[i,]+Ybar)
  plot(t,energy_data_day[i+1,],type='l',xlab='t',ylab='prediction of consumption',main=paste("i=",i),col='gray',lwd=2,ylim=range(energy_data_day))
  points(t,Yhattestrec,type='l',col='lightgreen',lwd=2)
}
legend('topleft',c("consumption","prediction"),lty=c(1,1),lwd=c(2,2),col=c('gray','lightgreen'))




source("../../PCA.R") 
source("../../functions_estim.R")

#################################################################################################################
# Extract of an initial file written by Dominik Liebl (importation of data)                                     #
#################################################################################################################
data.all <- read.csv(file="Data_Jan_2006_Sep_2008.csv", sep=";", header=TRUE)

## data-selector:=======================
peakload    <- rep(c(rep(FALSE,8),rep(TRUE,12),rep(FALSE,4)), (dim(data.all)[1]/24))
data.peakh  <- data.all[peakload,]
data.offph  <- data.all

## data.selector can be either data.selector="all", "peakh", or "offph"
data.list       <- switch("all",
                          all   = list("data"=data.all,  "T"=dim(data.all)[1]/24,   "N"=24),
                          peakh = list("data"=data.peakh,"T"=dim(data.peakh)[1]/12, "N"=12),
                          offph = list("data"=data.offph,"T"=dim(data.offph)[1]/24, "N"=24))
##=======================================                     

## day-index variable
first.day   <- min(data.list$data$DayIndex)
last.day    <- max(data.list$data$DayIndex)
days        <- c(first.day:last.day)

## Selection of non-workdays:
source("holidays_and_brueckentage.R")
Sys.getlocale("LC_TIME")
Sys.setlocale("LC_TIME", "C")
Da.vec <- matrix(data.all$Date, nrow=24)[1,]
Da.vec <- as.Date(Da.vec, "%d.%m.%Y")
feier  <- as.Date(feiertage, "%d.%m.%Y")

DUMMYSU    <- as.numeric(weekdays(Da.vec)=="Sunday")
DUMMYSA    <- as.numeric(weekdays(Da.vec)=="Saturday")
DUMMYMO    <- as.numeric(weekdays(Da.vec)=="Monday")

DUMMYFEIER    <- as.numeric(!is.na(match(Da.vec,feier)))
DUMMYPSEUDOMO <- as.numeric(!is.na(match(Da.vec,feier+1)))

non.workdays  <- c(1:length((Da.vec)))[c(as.logical(DUMMYSU)|as.logical(DUMMYSA)|as.logical(DUMMYFEIER))]

## Taking non-workdays from the sample
tmp  <- match(non.workdays, days)
tmp  <- tmp[!is.na(tmp)]
days <- days[-tmp]

## Number of days 'T' and hours 'N':
T <- length(days)
N <- data.list$N

## Container for data
RL.Prices <- matrix(NA, N, T)
RL.Wind   <- matrix(NA, N, T)


j <- 1
for(t in (days)){
  RL.Prices[,j]              <- data.list$data$Price[data.list$data$DayIndex==t]
  ##
  RL.Wind[,j]                <- data.list$data$Windinfeed[data.list$data$DayIndex==t]
  ##
  j <- j+1
}

#################################################################################################################
# Estimation of slope operator (written by G. Chagny, A. Meynaoui and A. Roche)                                 #
#################################################################################################################

# Remove outliers

maxs = apply(RL.Prices,2,max)
Q1 = quantile(maxs,0.25)
Q3 = quantile(maxs,0.75)
indext = which(maxs<=Q3+1.5*(Q3-Q1))

# Covariable X: eolien infeed

Xbar = rowMeans(RL.Wind[,indext])
X = t(RL.Wind[,indext]-Xbar)

# Variable Y to predict: electricity price

Ybar = rowMeans(log(RL.Prices[,indext]+1))
Y = t(log(RL.Prices[,indext]+1)-Ybar)

t = 1:24
n = nrow(X)
p = ncol(X)



par(mfrow=c(1,2))
plot(range(t),range(RL.Wind),type='n',xlab='time(hour)',ylab='wind infeed (kWh)',main="Original data")
for (i in indext){
  points(t,RL.Wind[,i],type='l',lty=2,col='gray')
}
points(t,Xbar,type='l',col=2,lty=2,lwd=2)
plot(range(t),range(X),type='n',xlab='time(hour)',ylab='',main='Centered data')
for (i in 1:n){
  points(t,X[i,],type='l',col='lightblue')
}
abline(h=0,col='red',lwd=2)



par(mfrow=c(1,2))
plot(range(t),range(RL.Prices[,indext]),type='n',xlab='time(hour)',ylab='price (€/kWh)',main="Original data")
for (i in indext){
  points(t,RL.Prices[,i],type='l',lty=2,col='gray')
}
points(t,rowMeans(RL.Prices[,indext]),type='l',col=2,lty=2,lwd=2)
plot(range(t),range(Y),type='n',xlab='time(hour)',ylab='',main='Centered data')
for (i in 1:n){
  points(t,Y[i,],type='l',col="lightblue")
}
abline(h=0,col='red',lwd=2)


# Cross-validation
kappa = 0.6
Nn = 15


res.CV <- CVpred(X,Y,kappa,Nn)

par(mfrow=c(1,1))
barplot(prop.table(table(res.CV$dim_selec)),col='blue',xlab="Dimension selected",ylab='Frequency')
hist(res.CV$pred_error,freq=F,col='orange',xlab='L^2 prediction error',main='')


par(mfrow=c(1,3))
ord = sort.int(res.CV$pred_error,index.return = T)$ix
show = ord[c(1,round(n/2),n)]
par(mfrow=c(1,3))
for (i in show){
  plot(t,Y[i,],type='l',xlab='t',ylab='Prediction of Y',col='lightblue',main=paste("i=",i),lwd=2)
  points(t,res.CV$Yhat[i,],type='l',col='darkgreen',lty=2,lwd=2)
}
legend('topleft',c(expression(Y[i]),expression(hat(Y)[i]^(-i))),lty=c(1,2),lwd=c(2,2),col=c('lightblue','darkgreen'))

# return to the electricity prices 

for (i in show){
  Yhattestrec = exp(res.CV$Yhat[i,]+Ybar)-1
  plot(t,RL.Prices[,indext[i]],type='l',xlab='t',ylab='prediction of consumption',main=paste("i=",i),col='gray',lwd=2,ylim=range(RL.Prices[,indext]))
  points(t,Yhattestrec,type='l',col='lightgreen',lwd=2)
}
legend('topleft',c("consumption","prediction"),lty=c(1,1),lwd=c(2,2),col=c('gray','lightgreen'))


library(latex2exp)
kappa <- c(10^seq(-25,-1,by = 1),seq(0.2,5,by=0.2)) #Valeurs de kappa = 8(1 + delta) suivant les notations de l'article
vec <- 26:35
kappa <- kappa[vec]

Error1 <- read.table("Pred-Error-Ex1.txt",sep = "",dec = ".")  
Error2 <- read.table("Pred-Error-Ex2.txt",sep = "",dec = ".") 
Error3 <- read.table("Pred-Error-Ex3.txt",sep = "",dec = ".")

Error1 <- Error1[vec,1]   
Error2 <- Error2[vec,1] 
Error3 <- Error3[vec,1]

Dim1 <- read.table("Optim-Dim-Ex1.txt",sep = "",dec = ".")
Dim2 <- read.table("Optim-Dim-Ex2.txt",sep = "",dec = ".")
Dim3 <- read.table("Optim-Dim-Ex3.txt",sep = "",dec = ".")

Dim1 <- Dim1[vec,1]   
Dim2 <- Dim2[vec,1]
Dim3 <- Dim3[vec,1]

plot(kappa,Error1*10^5,type = 'l',lwd=1.5,col='red', 
     main=latex2exp::TeX("Model $(i)$"),
     xlab = latex2exp::TeX("Value of $\\kappa$"), 
     ylab = latex2exp::TeX("$10^5 \\times$ Empirical MSPE"))

plot(kappa,Error2*10^3,type = 'l',lwd=1.5,col='red', 
     main=latex2exp::TeX("Model $(ii)$"),
     xlab = latex2exp::TeX("Value of $\\kappa$"), 
     ylab = latex2exp::TeX("$10^3 \\times$ Empirical MSPE"))

plot(kappa,Error3,type = 'l',lwd=1.5,col='red', 
     main=latex2exp::TeX("Model $(iii)$"),
     xlab = latex2exp::TeX("Value of $\\kappa$"), 
     ylab = latex2exp::TeX("Empirical MSPE"))


plot(kappa,Dim1,type = 'l',lwd=1.5,col='blue',
     main=latex2exp::TeX("Model $(i)$"),
     xlab = latex2exp::TeX("Value of $\\kappa$"),
     ylab = "Mean Optimal Dimension")
plot(kappa,Dim2,type = 'l',lwd=1.5,col='blue',
     main=latex2exp::TeX("Model $(ii)$"),
     xlab = latex2exp::TeX("Value of $\\kappa$"),
     ylab = "Mean Optimal Dimension")
plot(kappa,Dim3,type = 'l',lwd=1.5,col='blue',
     main=latex2exp::TeX("Model $(iii)$"),
     xlab = latex2exp::TeX("Value of $\\kappa$"),
     ylab = "Mean Optimal Dimension")


#Les tracÃ©es des Boxplots
Box1 <- read.table("Boxplot-Error1.txt",sep = "",dec = ".")
Box2 <- read.table("Boxplot-Error2.txt",sep = "",dec = ".")
Box3 <- read.table("Boxplot-Error3.txt",sep = "",dec = ".")

mean1 <- 1.52926808562334
mean2 <- 1.50979486356568
mean3 <- 1.14980630337086

boxplot(Box1[,1]*10^5,ylab = latex2exp::TeX("$10^5 \\times$ Empirical MSPE"),
        main=latex2exp::TeX("Model $(i)$"))
abline(h = mean1, col = "red",lwd = 2, lty = 2)
boxplot(Box2[,1]*10^3,ylab = latex2exp::TeX("$10^3 \\times$ Empirical MSPE"),
        main=latex2exp::TeX("Model $(ii)$"))
abline(h = mean2, col = "red",lwd = 2, lty = 2)
boxplot(Box3[,1],ylab = latex2exp::TeX("Empirical MSPE"),
        main=latex2exp::TeX("Model $(iii)$"))
abline(h = mean3, col = "red",lwd = 2, lty = 2)


#Tracés des boxplots en fonction de la taille 
Box1Sz <- read.table("Boxplot-Error1-Size.txt",sep = "",dec = ".")
Box2Sz <- read.table("Boxplot-Error2-Size.txt",sep = "",dec = ".")
Box3Sz <- read.table("Boxplot-Error3-Size.txt",sep = "",dec = ".")

Box1Sz <- as.matrix(Box1Sz) 
Box2Sz <- as.matrix(Box2Sz)
Box3Sz <- as.matrix(Box3Sz)

colnames(Box1Sz) <- as.character(c(200,400,600))
colnames(Box2Sz) <- as.character(c(200,400,600))
colnames(Box3Sz) <- as.character(c(200,400,600))

par(mar = c(5.1, 4.3, 4.1, 2.1))
xi <- 0.61
x0 <- xi + 1.0*0:(ncol(Box1Sz)-1) 
xf <- 1.39
x1 <- xf + 1.0*0:(ncol(Box1Sz)-1) 
y0 <- 10^5*apply(Box1Sz, 2, mean)

boxplot(Box1Sz*10^5,ylab = latex2exp::TeX("$10^5 \\times$ Empirical MSPE"),
        xlab = latex2exp::TeX("Sample size"), main=latex2exp::TeX("Model $(i)$"))
segments(x0 = x0, y0 = y0, x1 = x1, y1 = y0,
         col = "red", lwd = 2,lty = 5)

y0 <- 10^3*apply(Box2Sz, 2, mean)

boxplot(Box2Sz*10^3,ylab = latex2exp::TeX("$10^3 \\times$ Empirical MSPE"),
        xlab = latex2exp::TeX("Sample size"), main=latex2exp::TeX("Model $(ii)$"))
segments(x0 = x0, y0 = y0, x1 = x1, y1 = y0,
         col = "red", lwd = 2,lty = 5)

y0 <- apply(Box3Sz, 2, mean)

boxplot(Box3Sz,ylab = latex2exp::TeX("Empirical MSPE"),
        xlab = latex2exp::TeX("Sample size"), main=latex2exp::TeX("Model $(iii)$"))
segments(x0 = x0, y0 = y0, x1 = x1, y1 = y0,
         col = "red", lwd = 2,lty = 5)





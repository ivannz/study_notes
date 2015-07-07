setwd("~/Documents/Universite/Moscow/Teaching/Statistical Learning/R/")
rm(list=ls())

#########################################################
### Reference: 7.8 Lab: Non-linear Modeling in JWHT #####
#########################################################
library(ISLR)
attach(Wage)
dim(Wage)
str(Wage)
plot(Wage$age,Wage$wage) #Visualise the data


#1. POLYNOMIAL REGRESSION
fit=lm(wage~poly(age,4),data=Wage) 
coef(summary(fit))
#Compare with
fit1=lm(wage~age+I(age^2)+I(age^3)+I(age^4),data=Wage) 
coef(summary(fit1))

#The poly function returns orthogonal polynomials:
x<-seq(-5,5,by=.1)
orthpx<-poly(x,degree=5)
plot(x,orthpx[,1],type="l")
lines(x,orthpx[,2],col="magenta")
lines(x,orthpx[,3],col="blue")
lines(x,orthpx[,4],col="red")
lines(x,orthpx[,5],col="cyan")
#compare with raw polynomials
px<-poly(x,degree=7,raw=T)
plot(x,px[,1],type="l")
lines(x,px[,2])
lines(x,px[,3],col="blue")
lines(x,px[,4])
lines(x,px[,5],col="red")
lines(x,px[,7],col="cyan")

fit2=lm(wage~poly(age,4,raw=T),data=Wage) #Uses raw, not othogonal polynomials
coef(summary(fit2))

#Making predictions
agelims=range(age)
age.grid=seq(from=agelims[1],to=agelims[2])
preds=predict(fit,newdata=list(age=age.grid),se=TRUE)
se.bands=cbind(preds$fit+2*preds$se.fit,preds$fit-2*preds$se.fit)
par(mar=c(4.5,4.5,1,1), oma=c(0,0,4,0))
plot(age,wage,xlim=agelims ,cex=.5,col="darkgrey")
title("Degree-4 Polynomial", outer=T)
lines(age.grid,preds$fit,lwd=2,col="blue")
matlines(age.grid,se.bands,lwd=1,col="blue",lty=3)
#Selection of the order of the polynomial
fit.1=lm(wage~age,data=Wage)
fit.2=lm(wage~poly(age,2),data=Wage)
fit.3=lm(wage~poly(age,3),data=Wage)
fit.4=lm(wage~poly(age,4),data=Wage)
fit.5=lm(wage~poly(age,5),data=Wage)
anova(fit.1,fit.2,fit.3,fit.4,fit.5)
#Alternatively, use the fact that poly returns orthogomal polynomials
coef(summary(fit.5))  #p-values are the same

# Wild behaviour at the boundaries as the degree of the polynomial increases
fit15=lm(wage~poly(age,15),data=Wage)
preds15=predict(fit15,newdata=list(age=age.grid),se=TRUE)
lines(age.grid,preds15$fit,lwd=2,col="magenta")

#2. SPLINE
library(splines)
fit=lm(wage~bs(age,knots=c(25,40,60)),data=Wage)
pred=predict(fit,newdata=list(age=age.grid),se=T)
plot(age,wage,col="gray")
lines(age.grid,pred$fit,lwd=2)
lines(age.grid,pred$fit+2*pred$se ,lty="dashed")
lines(age.grid,pred$fit-2*pred$se ,lty="dashed")

dim(bs(age,knots=c(25,40,60)))
dim(bs(age,df=6)) #Use the df function to specify the number of deg of freedom
str(bs(age,df=6))
attr(bs(age,df=6),"knots")


#3. NATURAL SPLINE
fit2=lm(wage~ns(age,df=4),data=Wage)
pred2=predict(fit2,newdata=list(age=age.grid),se=T)
lines(age.grid, pred2$fit,col="red",lwd=2)


#4. SMOOTHING SPLINE
fit8=smooth.spline(age,wage,df=8)
fit16=smooth.spline(age,wage,df=16)
fit32=smooth.spline(age,wage,df=32)
plot(age,wage,xlim=agelims ,cex=.5,col="darkgrey")
title (" Smoothing  Spline ")
lines(fit8,col="lightblue",lwd=2)
lines(fit16,col="red",lwd=2)
lines(fit32,col="black",lwd=2)

library(sfsmisc)
x=age; y=wage 
S8<-hatMat(x, trace= FALSE, pred.sm = function(x, y, ...) predict(smooth.spline(x, y, df=8), x)$y)
S16<-hatMat(x, trace= FALSE, pred.sm = function(x, y, ...) predict(smooth.spline(x, y, df=16), x)$y)
eeS8=eigen(S8); eval8=eeS8$values; evec8=eeS8$vectors
eeS16=eigen(S16); eval16=eeS16$values; evec16=eeS16$vectors
# a. Plot of the eigenvalues
par(oma=c(0,0,3,0))
plot(eval8[1:30], pch=19, col="red", type="o", xlab="index", ylab="Eigenvalue")
lines(eval16[1:30], pch=19, col="blue", type="o", xlab="index", ylab="Eigenvalue of the smoother matrix S")
legend(22,1,c("df=8", "df=16"),pch=c(19,19), col=c("red","blue"))
title("Eigenvalue of the smoother matrix S")
# b. Plot of the eigenvectors
par(mfrow=c(2,3),mar=c(1,1,1,1))
plot(age,evec16[,2])
plot(age,evec16[,3])
plot(age,evec16[,4])
plot(age,evec16[,5])
plot(age,evec16[,8])
plot(age,evec16[,10])
dev.off()



#########################################################
############### Bias Variance Tradeoff ##################
#########################################################
# p.158 in THF
x<-runif(100)
x<-sort(x)
y<-sin(12*(x+0.2))/(x+0.2)+rnorm(100)
plot(x,y,pch=19, xlab="X", ylab="y", main="Smoothing Spline")
par(new=T)
plot(sort(x),sin(12*(x+0.2))/(x+0.2),axes=F,col="purple",type="l",lwd=2,ylab="",xlab="")

#df=5
fit5=smooth.spline(x,y,df=5)
lines(fit5, col="lightblue", lwd=3)
S5<-hatMat(x, trace= FALSE, pred.sm = function(x, y, ...) predict(smooth.spline(x, y, df=5), x)$y)
Sigma5=S5%*%S5; v5=diag(Sigma5)
lines(fit5$x,fit5$y+2*sqrt(v5) ,lty="dashed", col="lightblue")
lines(fit5$x,fit5$y-2*sqrt(v5) ,lty="dashed", col="lightblue")
#eeS5=eigen(S5); eval5=eeS5$values; evec5=eeS5$vectors
#plot(x,evec5[,3])

#df=9
fit9=smooth.spline(x,y,df=9)
plot(x,y,pch=19, xlab="X", ylab="y", main="Smoothing Spline")
par(new=T)
plot(sort(x),sin(12*(x+0.2))/(x+0.2),axes=F,col="purple",type="l",lwd=2,ylab="",xlab="")
lines(fit9,col="green",lwd=3)
S9<-hatMat(x, trace= FALSE, pred.sm = function(x, y, ...) predict(smooth.spline(x, y, df=9), x)$y)
Sigma9=S9%*%S9; v9=diag(Sigma9)
lines(fit9$x,fit9$y+2*sqrt(v9) ,lty="dashed", col="green")
lines(fit9$x,fit9$y-2*sqrt(v9) ,lty="dashed", col="green")

#df=15
fit15=smooth.spline(x,y,df=15)
plot(x,y,pch=19, xlab="X", ylab="y", main="Smoothing Spline")
par(new=T)
plot(sort(x),sin(12*(x+0.2))/(x+0.2),axes=F,col="purple",type="l",lwd=2,ylab="",xlab="")
lines(fit15,col="magenta",lwd=3)
S15<-hatMat(x, trace= FALSE, pred.sm = function(x, y, ...) predict(smooth.spline(x, y, df=15), x)$y)
Sigma15=S15%*%S15; v15=diag(Sigma15)
lines(fit15$x,fit15$y+2*sqrt(v15) ,lty="dashed", ,col="magenta")
lines(fit15$x,fit15$y-2*sqrt(v15) ,lty="dashed", ,col="magenta")



setwd("~/Documents/Universite/Moscow/Teaching/Statistical Learning/R/")
rm(list=ls())

##################################################################
# Reference: Lab 5.3 in JWHT: Cross-Validation and the Bootstrap #
##################################################################
# We use the Auto data set


##################################################################
# 1. Validation set approach
library(ISLR)
attach(Auto)
str(Auto)

plot(horsepower, mpg, pch=19)

# Select a random subset of 196 observations out of the original 392 observations
train=sample(392,196)

# Fit a linear model and compute the MSE on the hold out set
lm.fit=lm(mpg~horsepower, data=Auto, subset=train)
summary(lm.fit)
mean((mpg-predict(lm.fit, Auto))[-train]^2)

# Same procedure, but with a quadratic fit
lm.fit2=lm(mpg~poly(horsepower, 2), data=Auto, subset=train) 
mean((mpg-predict(lm.fit2, Auto))[-train]^2)
horsepowerlims=range(horsepower)
horsepower.grid=seq(from=horsepowerlims[1],to=horsepowerlims[2])
preds<- predict(lm.fit2, newdata=list(horsepower=horsepower.grid))
plot(horsepower, mpg, pch=19)
lines(horsepower.grid, preds, lwd=2, col="blue")

# Plot of the variability of the esti,ation of the test error using validation
val.error=rep(0,5)
train=sample(392,196)
for (i in 1:5){
  lm.fit=lm(mpg~poly(horsepower, i), data=Auto)
  val.error[i]= mean((mpg-predict(lm.fit, Auto))[-train]^2)
}
plot(1:5, val.error, type="l", ylim=c(15,26))
#lines(val.error)


##################################################################
# 2. LLOCV

# Use funcion glm instead of lm
glm.fit=glm(mpg~horsepower, data=Auto)
coef(glm.fit)
lm.fit=lm(mpg~horsepower, data=Auto) 
coef(lm.fit)

library(boot)
glm.fit=glm(mpg~horsepower ,data=Auto)
cv.err=cv.glm(Auto,glm.fit)
str(cv.err)
cv.err$delta

#Polynomial fit
cv.error=rep(0,5)
for (i in 1:5){
  glm.fit= glm(mpg~poly(horsepower, i), data=Auto)
  cv.error[i]= cv.glm(Auto, glm.fit)$delta[1]
}
plot(1:5, cv.error, pch=19)


##################################################################
# 3. K-fold cross validation

cv.error.10= rep(0,10)
for (i in 1:10){
  glm.fit=glm(mpg~poly(horsepower, i), data=Auto)
  cv.error.10[i]= cv.glm(Auto, glm.fit, K=10)$delta[1]
}
cv.error.10



##################################################################
# 5. Bootstrap
#
# We assess the accuracy of the coefficient estimates in simple LR in the Auto data
boot.fn= function(data, index)
return(coef(lm(mpg~horsepower, data=data, subset=index)))

boot.fn(Auto ,1:392) #Estimation using the full dataset

boot.fn(Auto, sample(392, 392, replace=TRUE))

boot(Auto, boot.fn, 1000)
summary(lm(mpg~horsepower, data=Auto))$coef

boot.fn=function(data,index)
coefficients(lm(mpg~horsepower+I(horsepower^2), data=data, subset=index))
boot(Auto, boot.fn, 1000)






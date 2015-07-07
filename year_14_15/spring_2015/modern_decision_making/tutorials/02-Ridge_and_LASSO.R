setwd("~/Documents/Universite/Moscow/Teaching/Statistical Learning/R/")
rm(list=ls())

# We use the credit data set
Credit<-read.csv("Credit.csv")
names(Credit)
str(Credit)

# There are 400 observations, for 12 variables. The `Balance` variables gives the 
# average credit card debt for a number of individuals. The predictors are both 
# quantitative, such as `Age`, `Cards` (number of credit cards), `Education` 
# (number of years of education), `Income` (in thousands of dollars), `Limit` 
# (credit limit), `Rating` (credit rating), and qualitative, `Gender`, `Student`, 
# `Status` (marital status) and `Ethnicity` (Caucasian, African American or Asian).

# A scatter plot for each pair of the quantitative variables can be obtained using 
# the command `pairs`
pairs(Balance~Age+Cards+Education+Income+Limit+Rating,Credit,col="brown")

#Collinearity might be an issue: it reduces the accuracy of the estimates.
plot(Credit$Age,Credit$Limit,pch=19)
lm.fit1=lm(Balance~Age+Limit,data=Credit)
summary(lm.fit1)
plot(Credit$Rating,Credit$Limit,pch=19)
lm.fit2=lm(Balance~Rating+Limit,data=Credit)
summary(lm.fit2)

library(car)
vif(lm.fit1)
vif(lm.fit2)
lm.fit3=lm(Balance~Age+Rating+Limit,data=Credit)
vif(lm.fit3)

# We can fit a LR model to predict the balance as a function of `Ethnicity` only.
lm.fit=lm(Balance~Ethnicity,data=Credit)
summary(lm.fit)

# When dealing with qualitative variables, R automatically creates dummy variables. 
# The baseline is African American with $531. It is estimated that the Asian category 
# will have $18 less debt than the African American category, and that the Caucasian 
# category will have $12 less debt than the African American category.

# Linear regression on all the predictors:
Credit<-Credit[,-1]  #Remove the first column
lm.fit0=lm(Balance~.,data=Credit)
summary(lm.fit0)

# To implement ridge regression and the lasso, we make use of the package `glmnet`.
library(glmnet)
  
# The `model.matrix` function automatically converts the qualitative variables into 
# numerical variables, needed with the `glmnet` function which can only take numerical 
# inputs. The [,-1] removes the first column of the matrix, corresponding to the intercept. 
# When using both ridge regression and the lasso, the intercept needs to be removed.  
x=model.matrix(Balance~.,Credit)[,-1]
y=Credit$Balance
head(x)
head(y)
grid=10^seq(5,-2,length=100)
ridge.fit=glmnet(x,y,alpha=0,lambda=grid)

# The argument `alpha` indicates which method is used. The value `alpha=0` corresponds 
# to ridge regression, and `alpha=1` to the lasso. In between values of `alpha` correspond 
# to an elastic-net penalty.

# By default, the function `glmnet` standardize the variables. Use `standardize=FALSE` 
# if you do not want glmnet to standardize.
args(glmnet)

# For each value of $\lambda$, `glmnet` returns a vector of estimated coefficients, 
# which can be accessed using `coef()`.
dim(coef(ridge.fit))
coef(ridge.fit)[,20]

#We expect the coefficients to be close to the LS coefficients for small $\lambda$. 
# And indeed,
coef(ridge.fit)[,100]

# As $\lambda$ increases, expect the size of the coefficients to shrink.
grid[20]
sqrt(sum(coef(ridge.fit)[-1,20]^2))
grid[60]
sqrt(sum(coef(ridge.fit)[-1,60]^2))

# We can plot the evolution of the value of the coefficients as a function of $\lambda$.
plot(grid,coef(ridge.fit)[2,],ylim = c(-8, 2), log="x", type="l",col="red", xlab="lambda (log scale)", ylab="ridge regression coefficient estimates")  #variable income
points(grid,coef(ridge.fit)[3,], log="x", type="l",col="black")  #variable limit
points(grid,coef(ridge.fit)[4,], log="x", type="l",col="magenta")  #variable rating
points(grid,coef(ridge.fit)[7,], log="x", type="l",col="blue")  #variable education
legend(100,-4,c("income","limit","rating","education"), lty=c(1,1,1,1), col=c("red","black","magenta","blue"))

# Alternatively, as $\lambda$ increases, the l-2 norm of the vector of ridge estimates decreases. We can thus plot the value of the coefficients as a function of $||\hat{\beta}^{ridge}_\lambda||_2/||\hat{\beta}^{ls}||_2$.
ratiol2=sqrt(colSums(coef(ridge.fit)[-1,]^2)/sum(coef(lm.fit0)[-1]^2))
plot(ratiol2,coef(ridge.fit)[2,],ylim = c(-8, 2), type="l",col="red", xlab="amount shrunk", ylab="ridge regression coefficient estimates")
points(ratiol2,coef(ridge.fit)[3,],type="l",col="black")  
points(ratiol2,coef(ridge.fit)[4,],type="l",col="magenta")  
points(ratiol2,coef(ridge.fit)[7,],type="l",col="blue")  
legend(0.1,-4,c("income","limit","rating","education"),lty=c(1,1,1,1), col=c("red","black","magenta","blue"))

# For the lasso, we proceed as before, except that we set `alpha=1`.
grid=10^seq(5,-3,length=100)
lasso.fit=glmnet(x,y,alpha=1,lambda=grid)

# For $\lambda=0$, we obtain the LS estimate:
coef(lasso.fit)[,100]

# We can plot the coefficient estimates as a function of 
# $||\hat{\beta}^{lasso}_\lambda||_1/||\hat{\beta}^{ls}||_1$ since the lasso uses an 
# l-1 penalty.
ratiol1=sqrt(colSums(abs(coef(ridge.fit)[-1,]))/sum(abs(coef(lm.fit0)[-1])))
plot(ratiol1,coef(lasso.fit)[2,],ylim = c(-8, 2), type="l",col="red", xlab="amount shrunk", ylab="lasso coefficient estimates")
points(ratiol1,coef(lasso.fit)[3,],type="l",col="black")  
points(ratiol1,coef(lasso.fit)[4,],type="l",col="magenta")  
points(ratiol1,coef(lasso.fit)[7,],type="l",col="blue")  
legend(0.1,-4,c("income","limit","rating","education"),lty=c(1,1,1,1), col=c("red","black","magenta","blue"))

# Note on the coefficient estimates returned by `glmnet`. 
# Compare the values of the estimates returned by the following three scenarios.

# First, use `glmnet` with standardized inputs
mx<-colMeans(x); mx<- matrix(rep(mx,dim(x)[1]),nrow=dim(x)[1], ncol=dim(x)[2], byrow=TRUE)
sdx<-sqrt(colMeans((x-mx)^2)); sdx<- matrix(rep(sdx,dim(x)[1]),nrow=dim(x)[1], ncol=dim(x)[2], byrow=TRUE)
ridge.fit0=glmnet((x-mx)/sdx,y,alpha=0,lambda=25,standardize=FALSE)
coef(ridge.fit0)

# Then, consider unit variance columns
ridge.fit2=glmnet(x/sdx,y,alpha=0,lambda=25,standardize=FALSE)
coef(ridge.fit2)

# Finally, consider non-standardized variables
ridge.fit1=glmnet(x,y,alpha=0,lambda=25)
coef(ridge.fit1)

# What do you observe? How do you explain this? 
mean(y)  


# EXERCISE ##############################################################################
# Fit RR + lasso on the Hitters data. 
# Goal: Predict a baseball player’s Salary on the basis of various statistics 
# associated with performance in the previous year.
library(ISLR)

names(Hitters)
dim(Hitters)
head(Hitters)
sum(is.na(Hitters$Salary))
Hitters=na.omit(Hitters)
dim(Hitters)
head(Hitters)

pairs(Salary~AtBat+Hits+HmRun+CHits+Errors,Hitters,col="brown")

x=model.matrix(Salary~.,Hitters)[,-1]
head(x)
y=Hitters$Salary
head(y)

#RR
grid=10^seq(10,-2,length=100)
ridge.fit=glmnet(x,y,alpha=0,lambda=grid)
dim(coef(ridge.fit))
#Lasso
lasso.fit=glmnet(x,y,alpha=1,lambda=grid)
dim(coef(lasso.fit))

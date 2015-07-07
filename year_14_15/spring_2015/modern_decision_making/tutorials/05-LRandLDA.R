setwd("~/Documents/Universite/Moscow/Teaching/Statistical Learning/R/")
rm(list=ls())
library(ISLR)

# percentage returns for the S&P 500 stock index over 1250 days, 
# from the beginning of 2001 until the end of 2005.

names(Smarket)
head(Smarket)
dim(Smarket)
attach(Smarket)

######################################################################
# Reference: Lab 4.6 in JWHT: Logistic Regression, LDA, QDA, and KNN #
######################################################################

# We fit a LR model abd then a LDA/QDA in order to predict Direction
# using Lag1 through Lag5 and Volume.

# 1. LOGISTIC REGRESSION
glm.fit=glm(Direction~Lag1+Lag2+Lag3+Lag4+Lag5+Volume , data=Smarket, family=binomial)
summary(glm.fit)
coef(glm.fit)

glm.probs=predict(glm.fit, type="response")
glm.probs[1:10]
contrasts(Direction)

glm.pred=rep("Down",1250)
glm.pred[glm.probs >.5]="Up"
table(glm.pred, Direction)
1-mean(glm.pred==Direction ) #Estimate of the training error rate

train=(Year<2005) #Our training sample
Smarket.2005= Smarket[!train, ]
dim(Smarket.2005)
Direction.2005=Direction[!train]

glm.fit=glm(Direction~Lag1+Lag2+Lag3+Lag4+Lag5+Volume, data=Smarket, family=binomial, subset=train)
glm.probs=predict(glm.fit, Smarket.2005, type="response")

glm.pred=rep("Down",252)
glm.pred[glm.probs >.5]="Up"
table(glm.pred,Direction.2005)
1-mean(glm.pred==Direction.2005) #Estimate of the test error rate

# We refit the model with less predictors
glm.fit=glm(Direction~Lag1+Lag2,data=Smarket, family=binomial, subset=train)
glm.probs=predict(glm.fit, Smarket.2005, type="response")
glm.pred=rep("Down",252)
glm.pred[glm.probs >.5]="Up"
table(glm.pred,Direction.2005)
1-mean(glm.pred==Direction.2005) #Estimate of the test error rate
106/(106+76)


# 2. LDA
library(MASS)
lda.fit=lda(Direction~Lag1+Lag2, data=Smarket, subset=train)
lda.fit
plot(lda.fit)

lda.pred=predict(lda.fit, Smarket.2005)
names(lda.pred)
lda.class=lda.pred$class
table(lda.class, Direction.2005)
1-mean(lda.class==Direction.2005)

# Remark: Applying a 50 % threshold to the posterior probabilities allows us to 
# recreate the predictions contained in lda.pred$class
sum(lda.pred$posterior[,1]>=.5)
sum(lda.pred$posterior[,1]<.5)

# The posterior probability output by the model corresponds to the probability 
# that the market will decrease:
lda.pred$posterior[1:20,1]
lda.class[1:20]


# 3. QDA
qda.fit=qda(Direction~Lag1+Lag2, data=Smarket, subset=train)
qda.fit

qda.class=predict(qda.fit, Smarket.2005)$class
table(qda.class, Direction.2005)
mean(qda.class==Direction.2005)



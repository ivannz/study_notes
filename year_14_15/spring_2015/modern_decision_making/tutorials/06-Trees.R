setwd("~/Documents/Universite/Moscow/Teaching/Statistical Learning/R/")
rm(list=ls())

##################################################################
######## Reference: 8.3 Lab: Decision Trees in JWHT ##############
##################################################################

# I. CLASSIFICATION TREES

library(tree)   #The tree library is used to construct classification and regression trees
library(ISLR)

attach(Carseats)
# Description: A simulated data set containing sales of child car seats at 400 different stores.
# A data frame with 400 observations on the following 11 variables.
#
# Sales: Unit sales (in thousands) at each location
# CompPrice: Price charged by competitor at each location
# Income: Community income level (in thousands of dollars)
# ShelveLoc: A factor with levels Bad, Good and Medium indicating 
#            the quality of the shelving location for the car seats at each site
# .../...

str(Carseats)
summary(Carseats$Sales)

# Sales is a continuous variable, and so we begin by recoding it as a binary variable. 
# We use the ifelse() function to create a variable, called High, which takes on a value 
# of Yes if the Sales variable exceeds 8, and takes on a value of No otherwise.
High=ifelse(Sales<=8, "No", "Yes")
summary(High)
Carseats = data.frame(Carseats, High)


# We now use the tree() function to fit a classification tree in order to predict High 
# using all variables but Sales.
tree.carseats = tree(High~.-Sales, Carseats)
summary(tree.carseats) #Deviance is -2 sum_{m,k} n_mk log(p_mk)

plot(tree.carseats)
text(tree.carseats) # Top of the tree: ShelveLoc:ac -> meaning?
levels(Carseats$ShelveLoc)

plot(tree.carseats)
text(tree.carseats, pretty=0) #pretty=0 instructs R to include the category names for any qualitative predictors

tree.carseats

# Estimate of the test error
set.seed(2)
train = sample(1:nrow(Carseats), 200)
Carseats.test = Carseats[-train ,]
High.test = High[-train]
tree.carseats = tree(High~.-Sales, Carseats, subset=train)
tree.pred = predict(tree.carseats, Carseats.test, type="class") #the argument type="class" instructs R to return the actual class prediction
table(tree.pred, High.test)

# Pruning of the tree

# The function cv.tree() performs cross-validation in order to determine the optimal 
# level of tree complexity; cost complexity pruning is used in order to select a sequence
# of trees for consideration.
set.seed(3)
cv.carseats = cv.tree(tree.carseats, FUN=prune.misclass) #K-fold cross-validation, K=10 by default
names(cv.carseats)
# The cv.tree() function reports 
#   * the number of terminal nodes of each tree considered (size) 
#   * the corresponding cross-validation error rate (dev) [total # of misclassifications during CV]
#   * the value of the cost-complexity parameter used (k)
cv.carseats
# The value -Inf returns the original tree

par(mfrow=c(1,2))
plot(cv.carseats$size, cv.carseats$dev, type="b")
plot(cv.carseats$k, cv.carseats$dev, type="b")
dev.off()
plot(cv.carseats)

# We now apply the prune.misclass() function in order to prune the tree to obtain 
# the nine-node tree. 
# The argument 'best' is an integer requesting the size (i.e. number of terminal nodes) 
# of a specific subtree in the cost-complexity sequence to be returned.
prune.carseats=prune.misclass(tree.carseats, best=9)
# Alternatively: prune.carseats=prune.misclass(tree.carseats, k=1.75)

plot(prune.carseats)
text(prune.carseats, pretty=0)

# How well does this pruned tree perform on the test data set?
tree.pred=predict(prune.carseats, Carseats.test, type="class")
table(tree.pred, High.test)

# What if we increase the value of best?
prune.carseats=prune.misclass(tree.carseats, best=15)
plot(prune.carseats)
text(prune.carseats, pretty=0)
tree.pred=predict(prune.carseats, Carseats.test, type="class")
table(tree.pred, High.test)


# II. REGRESSION TREES

library(MASS)
attach(Boston)
str(Boston)
# The Boston data frame has 506 rows and 14 columns
#
# medv: median value of owner-occupied homes in \$1000s
# rm: average number of rooms per dwelling
# nox: nitrogen oxides concentration (parts per 10 million)
# crim: per capita crime rate by town
# lstat: lower status of the population (percent)
# .../...
#
train = sample(1:nrow(Boston), nrow(Boston)/2)
tree.boston=tree(medv~., Boston, subset=train)
summary(tree.boston)
# deviance = sum of squared errors for the tree.
plot(tree.boston)
text(tree.boston ,pretty=0)

# Pruning the tree 

cv.boston=cv.tree(tree.boston)
plot(cv.boston$size, cv.boston$dev, type='b')
# the most complex tree is selected by cross-validation. If we still wish to prune the tree..
prune.boston=prune.tree(tree.boston, best=5)
plot(prune.boston)
text(prune.boston, pretty=0)

# We use the unpruned tree to make predictions on the test set:
yhat = predict(tree.boston, newdata=Boston[-train ,])
boston.test = Boston[-train, "medv"]
plot(yhat, boston.test)
abline(0,1)
mean((yhat-boston.test)^2)



# III. BAGGING AND RANDOM FORESTS

library(randomForest)

# First, we use 13 predictors (mtry=13) -> do bagging str(Boston)
bag.boston=randomForest(medv~., data=Boston, subset=train, mtry=13, importance =TRUE)
bag.boston

yhat.bag = predict(bag.boston, newdata=Boston[-train,])
plot(yhat.bag, boston.test)
abline(0,1)
mean((yhat.bag-boston.test)^2)

# We can change the number of trees grown by randomForest() using the ntree argument
bag.boston=randomForest(medv~.,data=Boston, subset=train, mtry=13, ntree=25)
yhat.bag = predict(bag.boston, newdata=Boston[-train,])
mean((yhat.bag-boston.test)^2)

# Growing a random forest proceeds in exactly the same way, except that we use a smaller 
# value of the mtry argument. By default, randomForest() uses p/3 variables when building a 
# random forest of regression trees, and √p variables when building a random forest of 
# classification trees. Here we use mtry = 6.
rf.boston=randomForest(medv~., data=Boston, subset=train, mtry=6, importance =TRUE)
yhat.rf = predict(rf.boston, newdata=Boston[-train ,])
mean((yhat.rf-boston.test)^2)

importance(rf.boston, type=2)
# IncNodePurity: decrease in node purity (RSS for reg and Gini for classif)

varImpPlot(rf.boston, type=2)


# IV. BOOSTING

library(gbm)
# The function gbm has many parameters:
# * distribution: REG: "gaussian" (for minimizing squared error), "laplace" (for minimizing absolute error)
#                 CLAS: "bernoulli" or "adaboost"
# * n.trees: how many trees to grow?
# * interaction.depth: the depth of each tree
boost.boston=gbm(medv~., data=Boston[train,], distribution="gaussian", n.trees=5000, interaction.depth=4)
summary(boost.boston)
yhat.boost=predict(boost.boston, newdata=Boston[-train,], n.trees=5000)
mean((yhat.boost-boston.test)^2)






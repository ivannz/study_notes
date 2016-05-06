## R session
library( ISLR )

## Biary classification task
fit1 <- glm( Direction~Lag1+Lag2+Lag3+Lag4+Lag5+Volume, data = Smarket, family = binomial, subset = NULL )

X <- model.matrix( Direction~Lag1+Lag2+Lag3+Lag4+Lag5+Volume, data = Smarket )
fit <- glm.fit( X, Smarket$Direction, family = binomial() )




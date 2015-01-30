rm( list = ls( all.names = TRUE ) ) ; invisible( gc() );

data <- read.csv( "http://www-bcf.usc.edu/~gareth/ISL/Credit.csv" )
data <- data[, -1]
str( data )

library( ggplot2 )

ggplot( data ) +
	geom_point( aes( x = Income, y = Rating, colour = Ethnicity ) )

ggplot( data ) + geom_point( aes( Limit, Rating, colour = Ethnicity ) ) 

fit <- lm( Balance ~ Limit + Rating, data )
summary( fit )

## Strong linear relationship means that their coefficnient are more
##  undetermined. The problem of multicollinearity.

## The VIF -- variance inflation factor the ration of the variance of
##  the full model to the variance of the mode with one predictor.
library( car )

summary( fit <- lm( Balance ~ Limit + Rating + Ethnicity, data ) )
vif( fit )

## Moving on to implementation of the ridgre regression.
library( glmnet )

x <- model.matrix( Balance ~ . - 1, data )
y <- data$Balance

grid <- 10 ^ seq( 5, -2, length = 100 )
fit <- glmnet( x, y, alpha = 0, lambda = grid )

plot( grid, coef( fit )[7,], type = 'lines', col = "red" )

fit <- glmnet( x, y, alpha = 1, lambda = grid )
plot( grid, coef( fit )[7,], type = "line", col = "red" )

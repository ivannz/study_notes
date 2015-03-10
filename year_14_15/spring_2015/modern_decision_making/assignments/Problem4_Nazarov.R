rm( list = ls( all.names = TRUE ) ) ; invisible( gc( ) )

setwd( "~/study_notes/year_14_15/spring_2015/modern_decision_making/assignments/" )
# Credit <- read.csv( "http://www-bcf.usc.edu/~gareth/ISL/Credit.csv" )[ , -1 ]
Credit <- read.csv( "./Credit.csv" )[ , -1 ]

X <- model.matrix( Rating ~ . - 1, Credit )
Y <- Credit$Rating

sigma_1 <- function( rho, p = 4 ) rho ^ abs( outer( 1:p, 1:p,  `-` ) )
sigma_2 <- function( rho, p = 4 ) rho ^ abs( outer( 1:p, 1:p, `!=` ) )

ridge <- function( X, T, sigma = diag( ncol( X ) ), lambda = 10 ) {
	XX <- crossprod( X ) ; XT <- crossprod( X, T )
	sigma_inv <- solve( sigma )
	sapply( lambda, function( l ) crossprod( solve( XX + l * sigma_inv ), XT ) )
}

bayes_ridge <- function( X, T, rho = 0, lambda = 10, sig = NULL )
	lapply( rho, function( r ) ridge( X, T, sig( r, ncol( X ) ), lambda ) )

rho <- seq( -.99, .99, by = 0.01 )
lambda <- 10 ^ seq( -2, 2, length = 100 )
br1 <- bayes_ridge( X, Y, rho, lambda, sig = sigma_1 )
br2 <- bayes_ridge( X, Y, rho, lambda, sig = sigma_2 )

## Study the first coefficent
study <- function( br, i, main = NULL ) {
	cf <- sapply( br, `[`, i, 1:ncol( br[[ 1 ]] ) )
	image( lambda, rho, abs( cf ), useRaster = FALSE, col = heat.colors( 35 ), main = main )
	return( cf )
}

op <- par( mfcol = c( 2, 2 ), mar = c( 0, 0, 0, 0 ),
		mar = c( 2, 2, 2, 2 ), cex.axis = .5  )
	res <- study( br2, 1, main = colnames( X )[ 1 ] )
	res <- study( br2, 2, main = colnames( X )[ 2 ] )
	res <- study( br2, 3, main = colnames( X )[ 3 ] )
	res <- study( br2, 4, main = colnames( X )[ 4 ] )
par( op )

plot( grid, coef( res )[2,], type = "l",
	ylim = c( -8, 2 ), log = "x", col = "red",
	xlab = "lambda (log scale)",
	ylab = "ridge regression coefficient estimates" )

lines( grid, coef( res )[3,], col = "black" )
lines( grid, coef( res )[4,], col = "magenta" )
lines( grid, coef( res )[7,], col = "blue" )

legend( "bottomleft", c( "income", "limit", "rating", "education" ),
	lty = c( 1, 1, 1, 1 ), col = c( "red", "black", "magenta", "blue" ) )


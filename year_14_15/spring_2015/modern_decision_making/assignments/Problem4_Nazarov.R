rm( list = ls( all.names = TRUE ) ) ; invisible( gc( ) )

setwd( "~/study_notes/year_14_15/spring_2015/modern_decision_making/assignments/" )
# Credit <- read.csv( "http://www-bcf.usc.edu/~gareth/ISL/Credit.csv" )[ , -1 ]
Credit <- read.csv( "./Credit.csv" )[ , -1 ]

sigma_1 <- function( rho, p = 4 ) rho ^ abs( outer( 1:p, 1:p,  `-` ) )
sigma_2 <- function( rho, p = 4 ) rho ^ abs( outer( 1:p, 1:p, `!=` ) )

ridge <- function( X, T, rho = 0, lambda = 10, sigma = NULL ) {
	sigma <- if( is.function( sigma ) ) sigma else function( r, p ) diag( p )
## Precompute the matrix-matrix products
	XX <- crossprod( X ) ; XT <- crossprod( X, T )
## For each $\rho$
	sim <- lapply( rho, function( r ) {
##  compute the inverse of $\Sigma$
		sigma_inv <- solve( sigma( r, ncol( X ) ) )
##  and the ridge estimator.
		sapply( lambda, function( l ) crossprod( solve( XX + l * sigma_inv ), XT ) )
	} )
	# sim <- lapply( rho, function( r ) ridge( X, T, sig( r, ncol( X ) ), lambda ) )
## Present the results in by-coefficient format
	if( ncol( X ) > 1 ) {
		structure( names = colnames( X ), 
			lapply( 1 : ncol( X ), function( i )
				aperm( sapply( sim, `[`, i, seq_along( lambda ) ), c( 2, 1 ) ) ) )
	} else {
		structure( names = colnames( X ), 
			list( aperm( simplify2array( sim ), c( 2, 1 ) ) ) )
	}
}

## Study the absolute value of a coefficent by plotting it as a heat map
plt <- function( res, i )
	image( rho, lambda, abs( res[[ i ]] ), col = heat.colors( 12 ),
		main = names( res )[ i ] )

evl <- function( res, j, i ) {
	data <- lapply( j, function( x ) res[[x]][i,] )

	plot( range( lambda ), range( unlist( data ) ), type = "n",
		log = "x", xlab = "lambda (log scale)",
		ylab = "Bayes ridge regression coefficient estimates" )

	colours <- topo.colors( n = length( j ) )
	invisible( lapply( seq_along( data ), function( x ) {
		lines( lambda, data[[x]], col = colours[ x ] )
	} ) )
	legend( "bottomleft", names( res )[ j ], ncol = 2,
		lty = c( 1, 1, 1, 1 ), col = colours )
}

## First center and scale the input and the response
X <- scale( model.matrix( Balance ~ . - 1, Credit ) )
Y <- scale( Credit$Balance )

rho <- seq( -.99, .99, by = 0.01 )
lambda <- 10 ^ seq( -2, 2, length = 100 )
br1 <- ridge( X, Y, rho, lambda, sigma = sigma_1 )
br2 <- ridge( X, Y, rho, lambda, sigma = sigma_2 )

plt( br1, 1 )

op <- par( mfcol = c( 2, 2 ), mar = c( 0, 0, 0, 0 ),
		mar = c( 2, 2, 2, 2 ), cex.axis = .5  )
	evl( br2, 0 + 1:3, 87 )
	evl( br2, 3 + 1:3, 87 )
	evl( br2, 6 + 1:3, 87 )
	evl( br2, 9 + 1:3, 87 )
par( op )

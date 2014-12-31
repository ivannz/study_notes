fn <- function( x, lambda = 3, gamma = 6, theta = 0.9, epsilon = 0.2 ) {
## Main hypothesized distribution
	f <- lambda * exp( - lambda * x )
	F <- ( 1 - exp( - lambda * x ) )

## Contaminating admixture
	u <- ( 0.1 <= x & x <= 0.4 ) / 0.3
	U <- ifelse( 0.1 <= x & x <= 0.4, ( x - 0.1 ) / 0.3, 0 )

## Resulting density
	g_k <- theta * f + ( 1 - theta ) * u
	G_k <- theta * F + ( 1 - theta ) * U

## tail is exp( -\gamma \abs{ x } )
	G_k + g_k / gamma - 1 / ( 1 + epsilon )
}

a <- uniroot( fn, c( 0, 1e04 ), tol = 0.0001 )$root









rm( list = ls( all.names = TRUE ) ) ; invisible( gc( ) )

sigma <- function( rho, p = 4 ) {
	pow <- abs( (1-p):(p-1) )
	mat <- matrix( 0, p, p )
	for( i in 1:p )
		mat[ i, ] <- pow[ 1:p-i ]
}

rho = .97
p = 4


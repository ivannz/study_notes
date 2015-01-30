## A simple simulation to test whether the explicite formula for the
##  probability is close to the truth.
theor_pr <- function( n, k ) {
	if( n < k ) return( 0 );
	m <- seq( 0, k-1, 1 )
	lnCmk <- lgamma( k+1 ) - lgamma( k-m+1 ) - lgamma( m+1 )
	max( 0, min( 1, sum(
## k! {n,k} k^{-n} -- the chance that a given collection
##  admits a partition into ordered k subsets.
		exp( lnCmk + n * log( 1 - m/k ) ) * ( -1 ) ^ m ) ) )
}

simul_pr <- function( n, N = 1000 ) {
	if( n < k ) return( 0 );
	smpl <- replicate(
		sample( k, n, replace = TRUE ),
		n = N, simplify = FALSE )
## Count the number of distinct values in each replication
	distinct_vals <- sapply( lapply( smpl, unique ), length )
## Return the probability estimate
	sum( distinct_vals == k ) / N
}

## Run the simulation for the k possible label values
##  for as many n as possible
k <- 75
n <- seq( k-1, 1000 )

theor_prob <- sapply( n, theor_pr, k=k )
theor_prob2 <- sapply( n, theor_pr2, k=k )
simul_prob <- sapply( n, simul_pr )

plot( n, simul_prob, col = "blue", type = "l", lwd = 1 )
lines( n, theor_prob, col = "red", lwd = 2 )
lines( n, theor_prob2, col = "black", lwd = 2 )
lines( n, cumsum( theor_prob2 ), col = "green", lwd = 2 )

theor_pr2 <- function( n, k ) {
	if( n < k ) return( 0 );
	m <- seq( 0, k-2, 1 )
	lnCmk <- lgamma( k+1 ) - lgamma( k-m ) - lgamma( m+1 )
	max( 0, min( 1, sum(
		exp( lnCmk + ( n - 1 ) * log( 1 - ( m + 1 ) / k ) ) * ( -1 ) ^ m ) ) )
}

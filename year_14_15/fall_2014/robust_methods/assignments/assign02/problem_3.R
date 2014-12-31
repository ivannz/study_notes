rm( list = ls( all.names = TRUE) ) ; invisible( gc( ) )

setwd( "/users/user/Desktop/studies 2014-2015/Robust methods/assign02/tex" )

Sweave( file = "../assign02.Stex" )


B <- pnorm( .5 )
C <- 8 * ( 1 - B ) / ( 3 * pi )

## Define functions
## The actual distribution function
F <- function( x ) suppressWarnings(
  C * pi / 2 + ifelse( x < 0, C * atan( x ),
    B - 1 + ifelse( x <= 4, pnorm( ( x - 2 ) / 4 ),
      B + C/2 * atan( ( x - 4 )^2 ) ) ) )

## The inverse of the distribution function
F_inv <- function( u ) suppressWarnings(
  ifelse( u < C * pi / 2, tan( u / C - pi / 2 ),
    ifelse( u <= C * pi / 2 + 2 * B - 1,
      4 * qnorm( u - C * pi / 2 + 1 - B ) + 2,
        4 + sqrt( tan( 2 * ( u - C * pi / 2 - 2 * B + 1 ) / C ) ) ) ) )

## This function implements the Pearson-Tukey robust
##   mean estimator with 5 quantiles.
pt.mean <- function( x ) {
## The selected quantiles
  q <- c( 1, 4, 8, 12, 15 ) / 16
## ... are weighted according to:
  w <- c( 1, 1, 2,  1,  1 ) / 6
  sum( w * sort( x )[ round( q * length( x ) ) ] )
}

## An implementation of the Hodges-Lehman estimator
hl.mean <- function( x ) {
## Compute all pairwise means
  w <- outer( x, x, function( a, b ) .5*( a + b ) )
## Discard a half
  median( w[ outer( seq_along( x ), seq_along( x ), `>=` ) ] )
}

## Carry out a simulation study
sim_study <- function( nobs, samples,
  seed = NULL, fun = c( mean, pt.mean, hl.mean ) ) {
## Run the simulation
  if( is.numeric( seed ) ) set.seed( seed )
  generated <- replicate( n = samples, {
    F_inv( runif( n = nobs ) )
  }, simplify = FALSE )
## Compute the supplied functions of each replication
  names( fun ) <- as.character( substitute( fun )[-1] )
  result <- do.call( data.frame, args =
    lapply( fun, function( fn ) sapply( generated, fn ) ) )
}

## Estimate risk
risk <- function( data, risk )
  mean( do.call( risk, args = list( data ) ) )

## Risk functions: by default hypothesis a normal distribution with mean 2
##  and standard deviation 4. 
## Square Loss
r1 <- function( theta, theta_0 = 2 )
  abs( theta - theta_0 ) ^ 2

## risk well function
r2 <- function( theta, theta_0 = 2 )
  pmin( 1, abs( theta - theta_0 ) )

## L^1 risk with indistinguishability region
r3 <- function( theta, theta_0 = 2 )
  ifelse( abs( theta - theta_0 ) >= .3, abs( theta - theta_0 ), 0 )

## 
ss <- lapply( c( 1, 2, 4, 8, 16, 32, 64 ) * 10, function( nobs ) {
  mean_est <- sim_study( nobs = nobs, samples = 50 )
  sapply( c( r1 = r1, r2 = r2, r3 = r3 ),
    function( r ) sapply( mean_est, risk, risk = r ) )
} )


qplot( pt_mean )
qplot( hl_mean )



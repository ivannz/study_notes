rm( list = ls( all.names = TRUE) ) ; invisible( gc( ) )

setwd( "/users/user/Desktop/studies 2014-2015/Robust methods/assign01/tex" )

Sweave( file = "../assign01.Stex" )

## Form the matrix of factors
set.seed( -1234L )
## Generate observation errors (noise)
eta <- rnorm( N <- 12 )

factors <-
  expand.grid( x1 = c(-1, 0, 1), x2 = c(1, 2) )

F <- model.matrix( ~ 1 + x1 + x2 + x1*x2, factors )

## Set the experiment plan
obs <- c( 2, 3, 1, 3, 1, 2 )
R <- diag( 1, nrow( factors ) )[
  rep( seq_along( obs ), obs ), ]

## The vector of true parameters
C <- c( 0.3, -0.5, -0.5, 0.05 )

## The data generation process
Y <- ( X <- R %*% F ) %*% C + eta

## OLS of the unconstrained model Y ~ F
beta <- solve( crossprod( X ) ) %*% crossprod( X, Y )
S2_hat <- as.vector( crossprod( Y - X %*% beta ) )
####
beta ; S2_hat / ( N - 4 )

## Test whether x1 is insignificant
##  -- Drop columns 2 and 4
X_c <- X[ , c( -2, -4 ), drop = FALSE ]
beta_c <- solve( crossprod( X_c ) ) %*% crossprod( X_c, Y )
S2_star <- as.vector( crossprod( Y - X_c %*% beta_c ) )
F_stat <- ( S2_star - S2_hat ) * ( N - ncol( X ) ) / ( S2_hat * ( ncol( X ) - ncol( X_c ) ) )
pv <- pf( F_stat, ncol( X ) - ncol( X_c ), N - ncol( X ) )
####
beta_c ; S2_star
S2_star / ( N - ncol( X_c ) )
F_stat ; pv

## Test whether x2 is insignificant
##  -- Dropping columns 3 and 4
X_c <- X[ , c( -3, -4 ), drop = FALSE ]
beta_c <- solve( crossprod( X_c ) ) %*% crossprod( X_c, Y )
S2_star <- as.vector( crossprod( Y - X_c %*% beta_c ) )
F_stat <- ( S2_star - S2_hat ) * ( N - ncol( X ) ) / ( S2_hat * ( ncol( X ) - ncol( X_c ) ) )
pv <- pf( F_stat, ncol( X ) - ncol( X_c ), N - ncol( X ) )
####
beta_c ; S2_star
S2_star / ( N - ncol( X_c ) )
F_stat ; pv

## Test whether x1:x2 is insignificant
##  -- enforce c_4 = 0 constraint by dropping the fourth
##  variable altogether
X_c <- X[ , -4, drop = FALSE ]
##  -- Estimate the constrained model
beta_c <- solve( crossprod( X_c ) ) %*% crossprod( X_c, Y )
##  -- Compute its squared sum of errors
S2_star <- as.vector( crossprod( Y - X_c %*% beta_c ) )
##  -- Compute the likelihood ratio statistic and its P value
F_stat <- ( S2_star - S2_hat ) * ( N - ncol( X ) ) / ( S2_hat * ( ncol( X ) - ncol( X_c ) ) )
pv <- pf( F_stat, ncol( X ) - ncol( X_c ), N - ncol( X ) )
####
beta_c ; S2_star
S2_star / ( N - ncol( X_c ) )
F_stat ; pv


## Numerically unstable
Q <- matrix( c( 0, 1, -1, 0,
				1, 0,  0, 0 ), ncol = 4 )
q <- c( 0, 0 )
Q <- matrix( c( 0, 1, -1, 0 ), ncol = 4 )
q <- c( 0 )

M <- solve( Q %*% solve( crossprod( X ) )
	%*% t( Q ) ) %*% ( Q %*% beta - q )
S2_diff <- as.vector( t( Q %*% beta - q ) %*% M )
beta_c <- beta - solve( crossprod( X ) ) %*%  t( Q ) %*% M

F_stat <- S2_diff * ( N - ncol( X ) ) / ( S2_hat * nrow( Q ) )
pf( F_stat, nrow( Q ), N - ncol( X ) )

S2_star <- S2_hat + S2_diff

### Optimal paln
factors <-
  expand.grid( x1 = c(-1, 0, 1), x2 = c(1, 2) )

F <- model.matrix( ~ 1 + x1 + x2 + x1*x2, factors )

## Set the experiment plan
obs <- c( 2, 3, 1, 3, 1, 2 )
R <- diag( 1, nrow( factors ) )[
  rep( seq_along( obs ), obs ), ]
ei <- eigen( crossprod( R %*% F ) )
sqrt( ei$values[ 1 ] / ei$values[ 4 ] )

## Set the experiment plan
obs <- c( 2, 2, 2, 2, 2, 2 )
R <- diag( 1, nrow( factors ) )[
  rep( seq_along( obs ), obs ), ]
ei <- eigen( crossprod( R %*% F ) )
sqrt( ei$values[ 1 ] / ei$values[ 4 ] )

## Set the experiment plan
obs <- c( 3, 0, 3, 3, 0, 3 )
R <- diag( 1, nrow( factors ) )[
  rep( seq_along( obs ), obs ), ]
ei <- eigen( crossprod( R %*% F ) )
sqrt( ei$values[ 1 ] / ei$values[ 4 ] )

## Set the experiment plan
obs <- c( 5, 0, 1, 1, 0, 5 )
R <- diag( 1, nrow( factors ) )[
  rep( seq_along( obs ), obs ), ]
ei <- eigen( crossprod( R %*% F ) )
sqrt( ei$values[ 1 ] / ei$values[ 4 ] )



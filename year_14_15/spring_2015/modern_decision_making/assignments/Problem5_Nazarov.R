rm( list = ls( all.names = TRUE ) ) ; invisible( gc( ) )
devAskNewPage( ask = TRUE )

if( !tryCatch( { require( "RColorBrewer", quietly = TRUE ) ; TRUE }, error = function ( c ) FALSE ) ){
	install.packages( "RColorBrewer" )
	require( "RColorBrewer" )
}

library( MASS )
library( splines )

## Create a grid of distances
dis_seq <- do.call( seq, args = c( c( -1, 1 ) + range( Boston$dis ), list( by = .1 ) ) )

fit <- lm( nox ~ poly( dis, 3 ), data = Boston )
nox_dis <- predict( fit, newdata = list( dis = dis_seq ), se = TRUE )

## PLot the fitted polynomial
plot( Boston$dis, Boston$nox, type = "p", pch = "+",
	cex = 1, col = "gray75", xlab = "dis", ylab = "nox",
	main = "Cubic polynomial" )
lines( dis_seq, nox_dis$fit, col = "red", lwd = 2  )

## Plot ~95% confidense interval
lines( dis_seq, nox_dis$fit + 1.96 * nox_dis$se, col = "black", lwd = 1 )
lines( dis_seq, nox_dis$fit - 1.96 * nox_dis$se, col = "black", lwd = 1 )

## Show the regression output.
summary( fit )
## The regression is significant and the fit is adequate overall.

## Pick some colours and 
degrees <- seq( 1, 12, by = 1 )
names( degrees ) <-  paste0( "p", degrees ) 
colours <- brewer.pal( length( degrees ), "Paired" )

## Create the master plot
plot( Boston$dis, Boston$nox, type = "p", pch = "+",
	cex = 1, col = "gray75", xlab = "dis", ylab = "nox",
	main = "Cubic polynomial" )

## Fit a bunch of polynomials
RSS <- mapply( function( d, col ) {
	fit <- lm( nox ~ poly( dis, degreee = d ), data = Boston )
	lines( dis_seq, predict( fit, newdata = list( dis = dis_seq ) ), col = col, lwd = 2  )

	sum( fit$residuals * fit$residuals )
}, degrees, colours, USE.NAMES = TRUE )

## Add a legend 
legend( "topright", lwd = 2, col = colours, legend = names( degrees ), cex = 0.8, ncol = 2 )

## As the number of parameters grows the fit of the regression curve imporves, however
##  with dinimishing marginal reduction is the sum of squares.
print( RSS )
print( 1 - tail( RSS, -1 ) / head( RSS, -1 ) )

################################################################################
##                              CROSS VALIDATION                              ##
################################################################################
## Using some environment fiddling define a generic kfold cross validator of 
##  R's linear models.
kfold <-function( model, data = environment( model ), k = 10, seed = NULL, ... ) {
## Fix the model matrix for the cross validation
	X <- model.matrix( model, data = data, ... )
## Create the groupping labels: split into k COMPLETE groups (of equal size)
	groups <- rep( seq_len( k ), rep( nrow( X ) %/% k, k ) )
## Create a folding schedule
	schedule <- if( k < nrow( X ) ) {
## For a geniune k-folding generate a random permutaion
		if( !is.null( seed ) ) set.seed( seed )
## And split it according to the group labels
		split( sample.int( length( groups ) ), f = groups )
## No muscial reference here.
	} else {
## For leave-one-out cross validation, all observation are removed one by one.
		seq_len( nrow( X ) )
## C'mon, baby, do the LOOCV-omotion!
	}
## Perform the k-fold cross validation
	mse <- sapply( schedule, function( i ) {
## Exclude the fold and fit on the remaining data
		coef <- lm.fit( X[ -i, , drop = FALSE ], y = Boston$nox[ -i ] )$coefficients
## Compute the mean squared error of the model prediction on the excluded fold
		mean( ( Boston$nox[ i ] - X[ i, , drop = FALSE ] %*% coef )^2 )
	} )
## Compute the k-fold mean prediction squared error
	mean( mse )
}

################################################################################
## Create a set of models: this is meta!
poly_models <- lapply( degrees, function( d )
## Model sElection is Too metA. This is the only way to make R not
##  do its lazy evaluation, by addressing the environmento of this
##  function.
	as.formula( paste0( "nox ~ poly( dis, ", d, " )" ) ) )

## Do the most precise MSE first
loocv_mse <- sapply( poly_models, kfold, k = nrow( Boston ), data = Boston )

## The greater the number of folds, the less variance in the estimated MSE,
##  but slower is the computation. For practical reasons it turns out to
##  be efficient to use 20-fold x-validation.
folds <- c( 5, 10, 15, 20 ) ; names( folds ) <- paste0( folds, "-fold" )
## perform the k-fold cross validation
kfold_mse <- lapply( folds, function( k ) replicate( n = 20, 
	sapply( poly_models, kfold, k = k, data = Boston ) ) )

################################################################################
# Compute the mean values and the standard deviations of the fold MSE
avg_kfold_mse <- lapply( kfold_mse, function( mse ) apply( mse, 1, mean ) )
 sd_kfold_mse <- lapply( kfold_mse, function( mse ) apply( mse, 1, sd ) )

## Find the degree of the best polynomial in each type of cross validation
opt_loocv <- structure( names = "loocv", which.min( loocv_mse ) )
opt_kfold <- structure( names = names( folds ),
	sapply( avg_kfold_mse, which.min, USE.NAMES = FALSE ) )

################################################################################
## Create a plot
plot( range( degrees ), range( c( loocv_mse, unlist( avg_kfold_mse ) ) ), type = "n",
	xlab = "Degrees", ylab = "MSE (log)", main = "Cross validation", log = "y" )

## Brew some adequate colours
colours <- structure( c( "black", brewer.pal( length( folds ), "Paired" ) ),
	names = c( "loocv", names( folds ) ) )

## Plot the k-fold cross validation results
for( i in seq_along( folds ) ) {
	lines( degrees, avg_kfold_mse[[ i ]], type = "l", lwd = 2, col = colours[ i + 1 ] )
	abline( v = degrees[ opt_kfold[ i ] ], lwd = 1, col = colours[ i + 1 ] )
	# lines( degrees, avg_fold_mse_20 + 1.96 * sd_fold_mse_20, type = "l", lwd = 1, col = colours[ 2 ] )
}
## Plot the loocv results atop the others
lines( degrees, loocv_mse, type = "l", lwd = 2, col = colours[ 1 ]  )
abline( v = degrees[ opt_loocv ], lwd = 1, col = colours[ 1 ] )

## Add a lagend
legend( "topleft", lwd = 2, col = colours, legend = names( colours ), cex = 0.8, ncol = 2 )

################################################################################


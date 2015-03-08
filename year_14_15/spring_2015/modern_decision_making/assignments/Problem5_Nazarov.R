################################################################################
##                          PREAMBLE AND DECLARATIONS                         ##
################################################################################
rm( list = ls( all.names = TRUE ) ) ; invisible( gc( ) )
devAskNewPage( ask = TRUE )

## Color brewing functionality is required
if( !tryCatch( { require( "RColorBrewer", quietly = TRUE ) ; TRUE }, error = function ( c ) FALSE ) ){
	install.packages( "RColorBrewer" )
	require( "RColorBrewer" )
}

## Using some environment fiddling define a generic kfold cross validator of 
##  R's linear models.
lm.kfold <-function( model, data = environment( model ), k = 10, seed = NULL, ... ) {
## Check the number of folds
	if( k != -1 && k <= 1 ) warning( sprintf(
		"Invalid number of folds %d. Running LOOCV instead.\n", k ) )
## Fix the model matrix for the cross validation
	X <- model.matrix( model, data = data, ... )
## Create a folding schedule
	schedule <- if( 1 < k && k < nrow( X ) ) {
## Create the groupping labels: split into k COMPLETE groups (of equal size)
		groups <- rep( seq_len( k ), rep( nrow( X ) %/% k, k ) )
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
		mean( ( Boston$nox[ i ] - X[ i, , drop = FALSE ] %*% coef ) ^ 2 )
	} )
## Compute the k-fold mean prediction squared error
	mean( mse )
}

## Perform crossvalidation of linear models in batch
lm.batch_xvalidate <- function( models, folds = integer( 0 ),
	replications = 20, ..., .verbose = TRUE ) {
## Fix the extra arguments
	args <- list( ... )

## Leave-one-out cross validation
	if( .verbose ) cat( "LOOCV...\n" )
	loocv_mse <- sapply( X = models, FUN = lm.kfold, k = -1, data = Boston )

## perform the k-fold cross validation
	if( .verbose && length( folds ) ) cat( "k-fold...\n" )
	kfold_mse <- lapply( folds, extra = args, FUN = function( k, extra ) {
		replicate( n = replications,
			do.call( sapply, args = c( list(
				X = models, FUN = lm.kfold, k = k ), extra ) ) )
		} )

# Compute the mean values and the standard deviations of the fold MSE
	if( .verbose ) cat( "processing...\n" )
	avg_kfold_mse <- lapply( kfold_mse, function( mse ) apply( mse, 1, mean ) )
	 sd_kfold_mse <- lapply( kfold_mse, function( mse ) apply( mse, 1, sd ) )

## Find the degree of the best polynomial in each type of cross validation
	opt_loocv <- structure( names = "loocv", which.min( loocv_mse ) )
	opt_kfold <- structure( names = names( folds ),
		sapply( avg_kfold_mse, which.min, USE.NAMES = FALSE ) )

## Return a cross-validation results object
	if( .verbose ) cat( "done...\n" )
	return( list(
		loocv = list( mse = loocv_mse, opt = opt_loocv, se = numeric( 1 ), raw = loocv_mse ),
		kfold = list( mse = avg_kfold_mse, opt = opt_kfold, se = sd_kfold_mse, raw = kfold_mse )
	) )
}

## Get the Boston data
library( MASS )
## Use splines
library( splines )

## Create a grid of distances
dis_seq <- do.call( seq, args = c( c( -1, 1 ) + range( Boston$dis ), list( by = .1 ) ) )

################################################################################
##                             A CUBIC POLYNOMIAL                             ##
################################################################################
fit <- lm( nox ~ poly( dis, 3 ), data = Boston )
nox_dis <- predict( fit, newdata = list( dis = dis_seq ), se = TRUE )

## Plot the fitted polynomial
plot( Boston$dis, Boston$nox, type = "p", pch = "+",
	cex = 1, col = "gray75", xlab = "dis", ylab = "nox",
	main = "Cubic polynomial" )

## Plot the fitted curve
lines( dis_seq, nox_dis$fit, col = "red", lwd = 2  )

## Plot ~95% confidence interval of the true nox
lines( dis_seq, nox_dis$fit + 1.96 * nox_dis$se, col = "black", lwd = 1 )
lines( dis_seq, nox_dis$fit - 1.96 * nox_dis$se, col = "black", lwd = 1 )

## Show the regression output.
smry <- summary( fit )
## The regression is significant and the fit is adequate overall.
{
	print( coef( smry ), digits = 3 )
	cat( sprintf( "Cubic polynomial\n R^2:\t%.3f\n RSS:\t%.3f\n",
		smry$r.squared, sum( smry$residuals ^ 2 ) ) )
}

################################################################################
## We are intereseted in model selection among polynomials of different degrees
degrees <- seq( 1, 12, by = 1 ) ; names( degrees ) <-  paste0( "P", degrees ) 

## Model sElection is Too metA.
## Create a set of polynomial models: this is meta!
poly_models <- lapply( degrees, function( d ) {
## This is the only way to make R not do its lazy evaluation in
##  the environment of this function.
		as.formula( paste0( "nox ~ poly( dis, ", d, " )" ) )
	} )

## Create the master plot
plot( Boston$dis, Boston$nox, type = "p", pch = "+",
	cex = .75, col = "gray75", xlab = "dis", ylab = "nox",
	main = "Cubic polynomial" )

## Pick some colours and fit a bunch of polynomials
colours <- brewer.pal( length( poly_models ), "Paired" )
RSS <- mapply( function( model, col ) {
## Fit the model and plot the fitted curve
	fit <- lm( model, data = Boston )
	lines( dis_seq, predict( fit, newdata =
		list( dis = dis_seq ) ), col = col, lwd = 2  )
## Return the RSS
	sum( fit$residuals * fit$residuals )
}, poly_models, colours, USE.NAMES = TRUE )

## Add a legend 
legend( "topright", lwd = 2, col = colours,
	legend = names( degrees ), cex = 0.8, ncol = 2 )

## As the number of parameters grows the fit of the regression curve
##  improves, however with dinimishing marginal reduction is the sum
##  of squares.
print( RSS )

################################################################################
## The greater the number of folds, the less variance in the estimated MSE,
##  but slower is the computation. For practical reasons it turns out to
##  be efficient to use 20-fold x-validation.
folds <- c( 5, 10, 15, 20 ) ; names( folds ) <- paste0( folds, "-fold" )
poly_xv <- lm.batch_xvalidate( poly_models, folds = folds, data = Boston )

## Create a plot
plot( range( degrees ), range( c( poly_xv$loocv$mse, unlist( poly_xv$kfold$mse ) ) ),
	type = "n", xlab = "Degrees", ylab = "MSE (log)", main = "Cross validation", log = "y" )

## Brew some adequate colours
colours <- structure( names = c( "loocv", names( folds ) ),
	c( "black", brewer.pal( length( folds ), "Paired" ) ) )
## Plot the k-fold cross validation results
for( i in seq_along( folds ) ) {
	lines( degrees, poly_xv$kfold$mse[[ i ]], type = "l", lwd = 1, col = colours[ i + 1 ] )
	abline( v = degrees[ poly_xv$kfold$opt[ i ] ], lwd = 1, col = colours[ i + 1 ] )
}

## Plot the loocv results atop the others
lines( degrees, poly_xv$loocv$mse, type = "l", lwd = 2, col = colours[ 1 ]  )
abline( v = degrees[ poly_xv$loocv$opt ], lwd = 1, col = colours[ 1 ] )

## Add a lagend
legend( "topleft", cex = 0.8, ncol = 2, lwd = 2,
	col = colours, legend = names( colours ) )

## Demonstrate the k-fold MSE-optimal models
{
	cat( "Cross validation selected the folowing models:\n" )
	cat( paste0( "\t", names( folds ), " \t", sapply( poly_xv$kfold$opt,
		function( i ) paste0( deparse( poly_models[[ i ]] ) ) ), "\n" ) )
	cat( paste0( "\tLOOCV  \t", deparse( poly_models[[ poly_xv$loocv$opt ]] ), "\n" ) )
}

################################################################################
##                                   SPLINES                                  ##
################################################################################
fit <- lm( nox ~ bs( dis, df = 4, degree = 3 ), data = Boston )

## Setting 4 degrees of freedom for a cubic spline imposes 2nd order smoothing
##  constraints. At each region we have a 3rd degree polynomial with 4 parameters.
##   * Continuity removes one degree of freedom in the chice of the parametrs at
##     each knot;
##   * First order smoothness does the same, and so does the continuity of the 
##     second derivative.
## The effective df = (K + 1) * 4 - K * 3 - 1 = K + 3. This means that there is
##   just one knot, and thus the domain is plit in halves. This regression spline
##   is a piecewise 2-order smooth curve with one knot.
## By default the bs() selects the knots using the empirical quantiles.
##  With 3 degrees of freedom this spline would coincide with a simple 3rd degree
##  polynomial.
nox_dis <- predict( fit, newdata = list( dis = dis_seq ), se = TRUE )

## Plot the fitted polynomial
plot( Boston$dis, Boston$nox, type = "p", pch = "+",
	cex = 1, col = "gray75", xlab = "dis", ylab = "nox",
	main = "Cubic spline" )

## Plot the fitted spline and ~95% confidence interval of the true nox
lines( dis_seq, nox_dis$fit, col = "red", lwd = 2  )
lines( dis_seq, nox_dis$fit + 1.96 * nox_dis$se, col = "black", lwd = 1 )
lines( dis_seq, nox_dis$fit - 1.96 * nox_dis$se, col = "black", lwd = 1 )

## Show the regression output.
smry <- summary( fit )
## The regression is significant and the fit is adequate overall.
{
	print( coef( smry ), digits = 3 )
	cat( sprintf( "Regression spline\n R^2:\t%.3f\n RSS:\t%.3f\n",
		smry$r.squared, sum( smry$residuals ^ 2 ) ) )
}

################################################################################
degrees_of_freedom <- 2 + seq_len( 14 )
names( degrees_of_freedom ) <- paste0( "df=", degrees_of_freedom )
spline_models <- lapply( degrees_of_freedom, function( df )
	as.formula( paste0( "nox ~ bs( dis, df = ", df, ", degree = 3 )" ) ) )

## Create the master plot
plot( Boston$dis, Boston$nox, type = "p", pch = "+",
	cex = .75, col = "gray75", xlab = "dis", ylab = "nox",
	main = "Regression spline" )

## Pick some colours and fit a bunch of polynomials
colours <- brewer.pal( length( spline_models ), "Paired" )
RSS <- mapply( function( model, col ) {
## Fit the model and plot the fitted curve
	fit <- lm( model, data = Boston )
	lines( dis_seq, predict( fit, newdata =
		list( dis = dis_seq ) ), col = col, lwd = 2  )
## Return the RSS
	sum( fit$residuals * fit$residuals )
}, spline_models, colours, USE.NAMES = TRUE )

## Add a legend 
legend( "topright", lwd = 2, col = colours,
	legend = names( degrees_of_freedom ), cex = 0.8, ncol = 2 )

## As the number of parameters grows the fit residla variance should decrease
##  strangely, this doesn't hold here in the case of the spline. Probably due
##  to the regularization.
print( RSS )

################################################################################
## We are interested in the fittness of the spline with varying number of free
##  parameters.
folds <- c( 5, 10, 20 ) ; names( folds ) <- paste0( folds, "-fold" )
spline_xv <- lm.batch_xvalidate( spline_models, folds = folds, data = Boston )

plot( range( degrees_of_freedom ), range( c( spline_xv$loocv$mse, unlist( spline_xv$kfold$mse ) ) ),
	type = "n", xlab = "Degrees", ylab = "MSE (log)", main = "BS Cross validation", log = "y" )

colours <- structure( c( "black", brewer.pal( length( folds ), "Paired" ) ),
	names = c( "loocv", names( folds ) ) )

for( i in seq_along( folds ) ) {
	lines( degrees_of_freedom, spline_xv$kfold$mse[[ i ]], type = "l", lwd = 2, col = colours[ i + 1 ] )
	abline( v = degrees_of_freedom[ spline_xv$kfold$opt[ i ] ], lwd = 1, col = colours[ i + 1 ] )
}

lines( degrees_of_freedom, spline_xv$loocv$mse, type = "l", lwd = 2, col = colours[ 1 ]  )
abline( v = degrees_of_freedom[ spline_xv$loocv$opt ], lwd = 1, col = colours[ 1 ] )

legend( "topleft", lwd = 2, col = colours, legend = names( colours ), cex = 0.8, ncol = 2 )

{
	cat( "Cross validation selected the folowing models:\n" )
	cat( paste0( "\t", names( folds ), " \t", sapply( xv$kfold$opt, function( i ) paste0( deparse( spline_models[[ i ]] ) ) ), "\n" ) )
	cat( paste0( "\tLOOCV  \t", deparse( spline_models[[ xv$loocv$opt ]] ), "\n" ) )
}


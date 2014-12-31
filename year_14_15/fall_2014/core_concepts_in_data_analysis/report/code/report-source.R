## 0. Preamble
rm( list = ls(all.names = TRUE ) ) ; invisible( gc() )

library( e1071 )
library( ggplot2 )
library( xtable )
library( gridExtra )

theme0 <- function( ... ) theme_bw() + theme_minimal(base_size = 18) +
	theme( legend.position = "none" )

setwd( "/users/user/Desktop/studies 2014-2015/MIRKIN/report/" )

#### Load the data
## Read the names of the columns and then the data itself
wdbc.info <- scan( "./data/wdbc.cols", quiet = TRUE, sep = "", comment.char = "#", fill = TRUE,
	what = list( names = character( ), classes = character( ), description = character( ) ) )
names( wdbc.info$description ) <- wdbc.info$names

wdbc_data <- read.table( "./data/wdbc.data", colClasses = wdbc.info$classes,
	sep = ",", header = FALSE, col.names = wdbc.info$names )

## Features are computed from a digitized image of a fine needle aspirate (FNA) of a breast mass.
##  They describe characteristics of the cell nuclei present in the image. A few of the images
##  can be found at http://www.cs.wisc.edu/~street/images/

## 
## Short aliases for frequently used data fragments
wdbc.m <- wdbc_data[c(2, 03:12)]
wdbc.s <- wdbc_data[c(2, 13:22)] ; names( wdbc.s ) <- names( wdbc.m )
wdbc.w <- wdbc_data[c(2, 23:32)] ; names( wdbc.w ) <- names( wdbc.m )

## wdbc.m <- cbind( wdbc[ 2 ], wdbc[ 03:12 ] / wdbc[ 13:22 ] )
## mask <- apply( sapply( wdbc.m[ -1 ], is.finite ), 1, all )
## wdbc <- wdbc.m[ mask, ]
wdbc <- wdbc.w ##structure( wdbc.w, names = names( wdbc.m ) )


################################################################################
##  HW1
##  1. Each one or two or three person(s) find(s) a dataset of their liking on the internet
##  2. Starting writing a report. Title page, Explanation of the choice of the dataset.
##     Information of the dataset: features, number of entities, source address, examples of problems

## REPORT 

##  HW2
##  4. Compute the three centrality characteristics (mean, median, midrange), three spread
##     characteristics (variance, std, halfrange), mode and 25% percentile(s).

midrange <- function( x ) .5 * ( max( x ) + min( x ) )
halfrange <- function( x ) .5 * ( max( x ) - min( x ) )
histmode <- function( x, n = 30 ) {
## Split the data into the required number of bind of equal width
	HD <- hist( x, plot = FALSE, n = n )
## The mode is undefined for multimodal distributions
	HD$mids[ which.max( HD$density ) ]
}

## Generate a summary table
summary <- data.frame( stringsAsFactors = FALSE, row.names = names( wdbc[-1] )
	, mean = sapply( wdbc[-1], mean ), mode = sapply( wdbc[-1], histmode ), midrande = sapply( wdbc[-1], midrange )
	, variance = sapply( wdbc[-1], var ), std = sapply( wdbc[-1], sd ), halfrange = sapply( wdbc[-1], halfrange )
	, prc25 = sapply( wdbc[-1], quantile, probs = .25, names = FALSE )
		, median = sapply( wdbc[-1], median )
	, prc75 = sapply( wdbc[-1], quantile, probs = .75, names = FALSE ) )

## Print it as a properly fomatted TeX table. For reference see
##  http://bcb.dfci.harvard.edu/~aedin/courses/ReproducibleResearch/RR/exampleSweave.Rnw
print( xtable( summary, digits = 4, caption = "" ) )

##  3. Take one of the features in your dataset and display its histograms with different
##     numbers of bins so that their shapes differ from each other. Make a corresponding comment.
xlab <- function( x ) eval( parse( text = paste0( "expression( Delta == ", x, " ) " ) ) )
hist_plot <- ggplot( wdbc, aes( x = perimeter, fill = diagnosis ) ) + theme0( )
H <- lapply( c( 5, 15, 30 ), function( binwidth ) hist_plot + labs( x = xlab( binwidth ) ) +
	geom_histogram( aes( y = ..density.. ), colour = "black", binwidth = binwidth ) )

# print( do.call( grid.arrange, args = c( H, list( ncol = 3, nrow = 1 ) ) ) )
print( grid.arrange( H[[ 1 ]], H[[ 2 ]], H[[ 3 ]], ncol = 3, nrow = 1
	, main = "Histograms of cell perimeter for various width of bins" ) )

##  HW3
##  5. Find two features with more or less “linear-like” scatterplot.
##  6. Display their scatter-plot.
##  7. Build a linear regression of one of them over the other. Make a comment on the meaning of the slope.
##  8. Find the correlation and determinacy coefficients, and comment on the meaning of the latter.

ggplot( wdbc, aes( x = compactness, y = symmetry ) ) + theme0( ) + geom_point( aes( colour = diagnosis ) ) + geom_rug(size=0.1)

## Compute a simple bivariate linear regression
rho <- cor( wdbc$compactness, wdbc$symmetry )
beta <- rho * sd( wdbc$symmetry ) / sd( wdbc$compactness )
alpha <- mean( wdbc$symmetry ) - mean( wdbc$compactness ) * beta

## The R-squared (coefficient of determination) is in this particular case just the square of the correlation
Rsq <- rho^2

## Print the OLS estimation results
fit <- lm( symmetry ~ compactness, wdbc )
xtable( summary( fit ) )

##  9. Estimate the 95% confidence intervals of the slope, intercept, and correlation coefficient by using a 5000 strong bootstrap resampling technique. Make a comment.
bootstrap_denorm <- function( x, y ) {
## Draw observations with replacement
	bi <- sample( seq_along( y ), replace = TRUE )
## Get the sample correlation coefficient
	rho <- cor( bx <- x[ bi ], by <- y[ bi ] )
## Recover the slope and the intercept
	beta <- rho * sd( by ) / sd( bx )
	alpha <- mean( by ) - mean( bx ) * beta
## No need for standartization by s.e. just yet
	c( alpha = alpha, beta = beta, rho = rho )
}

## Run pairs bootstrap to estimate a bootstrap distribution of the parameters
result <- as.data.frame( do.call( rbind, replicate( n = 5000, bootstrap_denorm( y = wdbc$symmetry, x = wdbc$compactness ), simplify = FALSE ) ) )

## Display the results
xtable( do.call( rbind, lapply( result, quantile, probs = c( .025, .975 ) ) ) )

##  HW4
##  10. Take/make three nominal features at your dataset, x1, x2, x3. At least one of them should be made by categorizing a quantitative feature by using minima in its histogram.
##  11. Meaningfully take one of them, say x1, as the target and the other two as predictors. Build contingency tables (x2,x1) and (x3,x1).
##  12. Compute tables of conditional probability and Quetelet indexes. Comment of the large values.
##  13. Compute chi-squared values for the tables and comment of their meaning and of which of x2, x3 better correlates with x1.

## A function for categorizing data according to the supplied vector of split-points
categorize <- function( data, splits, digits = 3 ) {
## Arrange data into bins according to breaks
	bins <- .bincode( data, c( -Inf, splits, +Inf ) )

## Create suitable labels
	chsplits <- formatC( splits, digits = digits, width = 1L )
	labels <- c( paste0( "≤", head( chsplits, n = 1 ) ),
		paste0( "(", head( chsplits, n = -1 ), ", ", tail( chsplits, n = -1 ), "]" ),
		paste0( tail( chsplits, n = 1 ), "<" ) )

	factor( bins, seq_along( labels ), labels )
}

## Split into bins: split points are selected near local minima of the histogram
##  use xx <- hist( data, plot = FALSE, n = 40 )
# wdbc$cpt_bins <- categorize( wdbc$compactness, c( 0.2605, 0.2655, 0.2715 ) )
wdbc$cpt_bins <- categorize( wdbc$compactness, c( 0.11, 0.27, 0.45 ) )
wdbc$cvt_bins <- categorize( wdbc$concavity, c( 0.225, 0.475 ) )
wdbc$tex_bins <- categorize( wdbc$texture, c( 15, 22 ) )

## Construct a contingency table
ctable <- function( row, col, joint = FALSE, border = TRUE ) {
	core <- tapply( row, list( row, col ), length )
	core[ is.na( core ) ] <- 0
	ctab <- if( !border ) core else {
		rows <- cbind( core, sum = rowSums( core, na.rm = TRUE ) )
		rbind( rows, sum = colSums( rows, na.rm = TRUE ) )
	}
	if( joint ) ctab / sum( core ) else ctab
}

## Construct a table of Quetelet indices
quetelet <- function( row, col ) {
	joint <- ctable( row, col, joint = TRUE, border = TRUE )
## Compute the index using this formula Q = p(Gl & Hk) / ( p(Hk) p(Gl) ) - 1
	quetelet <- joint / outer( joint[ , ncol( joint ) ], joint[ nrow( joint ), ] ) - 1
## Remove the border
	quetelet[ -nrow( quetelet ), -ncol( quetelet ) ]
}

## Diagnosis ~ concavity contingency table and ChiSq test
ct_dia_cvt <- ctable( wdbc$diagnosis, wdbc$cvt_bins )
xtable( ct_dia_cvt, digits = 0 )
cs_dia_cvt <- chisq.test( ctable( wdbc$diagnosis, wdbc$cvt_bins, border = FALSE ) )

## Diagnosis ~ compactness contingency table and ChiSq test
ct_dia_cpt <- ctable( wdbc$diagnosis, wdbc$cpt_bins )
xtable( ct_dia_cpt, digits = 0 )
cs_dia_cpt <- chisq.test( ctable( wdbc$diagnosis, wdbc$cpt_bins, border = FALSE ) )

## Diagnosis ~ texture contingency table and ChiSq test
ct_dia_tex <- ctable( wdbc$diagnosis, wdbc$tex_bins )
xtable( ct_dia_tex, digits = 0 )
cs_dia_tex <- chisq.test( ctable( wdbc$diagnosis, wdbc$tex_bins, border = FALSE ) )

## Diagnosis ~ concavity Quetelet indices
qi_dia_cvt <- quetelet( wdbc$diagnosis, wdbc$cvt_bins )
xtable( qi_dia_cvt, digits = 3 )

## Diagnosis ~ compactness Quetelet indices
qi_dia_cpt <- quetelet( wdbc$diagnosis, wdbc$cpt_bins )
xtable( qi_dia_cpt, digits = 3 )

## Diagnosis ~ texture Quetelet indices
qi_dia_tex <- quetelet( wdbc$diagnosis, wdbc$tex_bins )
xtable( qi_dia_tex, digits = 3 )

ctable( wdbc$tex_bins, wdbc$cpt_bins )
quetelet( wdbc$tex_bins, wdbc$cpt_bins )
chisq.test( ctable( wdbc$tex_bins, wdbc$cpt_bins, border = FALSE ) )

## Describe transformation of the quetelet index
## Describe the Chi squared test and its relation to Quetelet index

##	HW5
##	Nonlinear regression (c.f. nonlinear_regression.docx)

## Regress 


##  HW 6
##  Multivariate regression using the nitty-gritty approach of projections



##	HW7
##	Adapt this code (nnbpm.m) to your dataset. Output results for different learning rates and different numbers of hidden nodes.



## Re-scale the data
resc <- function( X, a = -1, b = 1 ) {
## Force X into a matrix if needed
	if( !is.matrix( X ) )
		X <- as.matrix( X )
## Get the range of X
	xr <- apply( X, 2, range )
	xl <- xr[2, ] - xr[1, ]
## Normalize each column of X to [0,1]
	xu <- sweep( sweep( X, 2, xr[1,], `-` ), 2, xl, `/` )

## Re-scale
	xu * ( b - a ) + a
}

init_snn <- function( X, Y, H = 10 ) {
## Normalize the data 
	X <- cbind( resc( X, -1, 1 ), ..bias.. = 1 )
	Y <- resc( Y, -1, 1)

## Initialize the input-hidden layer weights
	W1 <- matrix( rnorm( ncol( X ) * H ), ncol( X ), H )
## ... and the hidden-output layer weights
	W2 <- matrix( rnorm( H * ncol( Y ) ), H, ncol( Y ) )

	return( list(
		data = list(
			input = X,
			output = Y ),
		model = list(
			layers = c( H ),
			weights = list( W1 = W1, W2 = W2 ) ) ) )
}

fit_snn <- function( object, mu = 0.001, maxiter = 10000, abs.error = 0.01 ) {
	W1 <- object$model$weights$W1
	W2 <- object$model$weights$W2
	X <- object$data$input
	Y <- object$data$output
## The main forward-backward pass loop
	niter <- 0
	repeat {
		sq <- sample( seq_len( nrow( X ) ) )
## Initialize mean absolute error to zero
		abs_err <- numeric( nrow( X ) )
		for( i in sq ) {
			x <- X[i,, drop = FALSE] ; y <- Y[i,, drop = FALSE]
## Forward pass
			OL1 <- 2 / ( 1 + exp(- x %*% W1) ) - 1
			OL2 <- OL1 %*% W2
## Accumulate the error
			ERR <- y - OL2
			abs_err[ i ] <- abs( ERR )

## Error back-prapagation.
## dW2 the gradient of the second layer
			dW2 <- -t( OL1 ) %*% ERR
			t1 <- W2 %*% t( ERR )
			t2 <- ( 1 - OL1 ) * ( 1 + OL1 ) / 2
			t3 <- t2 * t( t1 )
			dW1 <- -t( x ) %*% t3
## Update the weights of both layers
			W2 <- W2 - mu * dW2
			W1 <- W1 - mu * dW1
		}
## Stop if the required accuracy has been achieved
		if( mean( abs_err ) < abs.error ) break
## Or if the limit on the number of iterations has been exceeded
		if( niter > maxiter ) break
		niter <- niter + 1
	}

## Return a proper object
	new_model <- object$model
	new_model$weights <- list( W1 = W1, W2 = W2 )
	return( list(
		data = object$data,
		model = new_model,
		error = mean( abs_err ),
		niter = niter ) )
}

run_snn <- function( object, input ) {
	W1 <- object$model$weights$W1
	W2 <- object$model$weights$W2
	as.matrix( apply( input, 1, function( x ) {
		OL1 <- 2 / ( 1 + exp( -x %*% W1) ) - 1
		OL2 <- OL1 %*% W2
	} ) )
}

# X <- matrix( rnorm( 10*569 ), 569, 10 )
# Y <- rbinom( 569, 1, .5 ) 
X <- matrix( rnorm( 10*569 ), 569, 10 )
Y <- 

## Initialize the layered neural network model
mdl <- ( snn <- init_snn( X, Y, 10 ) )

## Split the sample in two: one for learning, the other for testing.
s0 <- sample( seq_len( nrow( X ) ), nrow( X ) * 2/3 )

## Restrict to the learning sample only
snn$data$input <- snn$data$input[s0,,drop=FALSE]
snn$data$output <- snn$data$output[s0,,drop=FALSE]

## Teach the neural network
snn <- fit_snn( snn, maxiter = 2000 )

## Test the neural network on the sample kept for cross-validation
tX <- mdl$data$input[-s0,,drop=FALSE]
tY <- mdl$data$output[-s0,,drop=FALSE]
oY <- run_snn( snn, tX )

  predicted <-
    factor( c( "B", "M" )[resc( oY >= 0, 1, 2)] )

  actual <-
    wdbc$diagnosis[ - train_sample ] #$

  library( tables )
  latex( tabular( actual ~ predicted), label = "tab:11",
    caption = "The confusion matrix of the neural network on the test sample" )

####################################################################################################
####################################################################################################
####################################################################################################
qplot( 100 * ( log( wdbc$compactness ) + log( 2 ) + .5 * log( pi ) ) )+theme0()
## aes( fill = ..density.. ) -- fill according to the values of density
## aes( fill = wdbc$dgn ) -- split-fill according to the count of binary variable

ggplot( wdbc, aes( x = texture, fill = diagnosis ) ) + geom_histogram( colour = "black" )
ggplot( wdbc, aes( x = compactness, fill = diagnosis ) ) + geom_histogram( colour = "black" )
ggplot( wdbc, aes( x = compactness, fill = diagnosis ) ) + geom_histogram( colour = "black" )
ggplot( wdbc, aes( x = smoothness, fill = diagnosis ) ) + geom_histogram( colour = "black" )

ggplot( wdbc, aes( x = texture, y = fractal_dimension ) ) + theme0( ) + geom_point( aes( colour = diagnosis ) )
ggplot( wdbc, aes( x = fractal_dimension, y = compactness ) ) + theme0( ) + geom_point( aes( colour = diagnosis ) )
ggplot( wdbc, aes( x = texture, y = compactness ) ) + theme0( ) + geom_point( aes( colour = diagnosis ) )
ggplot( wdbc, aes( x = perimeter, y = area ) ) + theme0( ) + geom_point( aes( colour = diagnosis ) )


invisible( lapply( names(wdbc[-1]), function( x ) {
	m <- ggplot( wdbc, eval( parse( text = paste0( "aes( x = ", x, " ) " ) ) ) ) +
		theme0() + geom_histogram( aes( fill = diagnosis ), colour = "black" )
	ggsave( file.path( ".", "data", "images", paste0( x, ".png" ) ), m )
} ) )


ggplot( wdbc, aes( x = sqrt( area ) / perimeter, fill = diagnosis ) ) + geom_histogram( colour = "black" )

m <- ggplot( wdbc, aes( x = compactness, fill = diagnosis ) ) +
	ggtitle( paste0( "Histogram of ", wdbc.info$description[ "compactness" ] ) ) +
	labs( x = wdbc.info$description[ "compactness" ] )
m + geom_histogram( colour = "black", binwidth = .05 )

ggplot( wdbc, aes( compactness, symmetry, fill = ..density.. ) ) + theme0( ) + stat_binhex( bins = 30 )

ggplot( result, aes( x = alpha ) ) + geom_histogram( ) + geom_vline( xintercept = quantile( alpha, probs = c( 0.025, 0.975) ) )
ggplot( result, aes( x = alpha ) ) + geom_histogram( aes( fill = as.factor( .bincode( alpha, c( -Inf, quantile( alpha, probs = c( 0.025, 0.975 ) ), +Inf ) ) ) ) )

## Better results are achieved with normalization by the s.e.
bootstrap_norm <- function( formula, data = NULL, N = 100, H = 0 ) {
## Setup a data frame for bootstrapping linear regression
	data <- model.frame( formula = formula, data = data )
## Replicate the bootrstap sampling the given number of times
	replicate( n = N, {
## Draw a random sample with replacement
		bi <- sample( nrow( data ), replace = TRUE )
## Estimate the model using lm
		bfit <- lm( formula, data = data[ bi, , drop = FALSE ] )
## Compute the bootstrapped t-statistic
		( coef( bfit ) - H ) / sqrt( diag( vcov( bfit ) ) )
	} )
}

zz <- bstrap( symmetry ~ compactness, data = wdbc, N = 100 )



# % In order to override the default behaviour of the plot aesthetics definitions, let's introduce a function \emph{xlab}, which returns an R expression object used for proper labelling of the x-axis of any histogram.
# % \begin{Scode}{echo=TRUE}
# % xlab <- function( x ) eval( parse( text =
# %   paste0( "expression( Delta == ", x, " ) " ) ) )
# % \end{Scode}


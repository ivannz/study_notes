# library( maps )
library( rgdal )
library( foreign )


setwd( "/users/user/Desktop/studies 2014-2015/MIRKIN/datasets/" )

ogrListLayers( "./buenosaires/buenosaires.shp" )
buenosaires <- readOGR(
	  "./buenosaires/buenosaires.shp"
	, "buenosaires", verbose = TRUE )

plot( buenosaires )

invisible( lapply( buenosaires@polygons, function( poly ) {
	polygon( poly@Polygons[[ 1 ]]@coords, col = "red", border = NA )
} ) )

xdbf<- read.dbf("./buenosaires/buenosaires.dbf")

################################################################################

id <- unique( replicate( n = 20, paste0( sample( letters, 2 ), collapse = "") ) )
data <- data.frame( cat = sample( id, 1000, replace = TRUE ), stringsAsFactors = FALSE )

## Convert categorical data to wide format using double indexed tapply()
.to_wide <- function( X, prefix = "" )
	tapply( rep( TRUE, length( X ) ), list( 1:length( X ), paste0( prefix, X ) ), c )

to_wide <- function(X, prefix = "" )
	{ W <- .to_wide( X, prefix ) ; W[ is.na( W ) ] <- FALSE ; W }

setwd( "/users/user/Desktop/studies 2014-2015/MIRKIN/datasets/" )
crx.raw <- read.table( "./credit scoring/crx.data", sep=",", na.strings = "?" )
names( crx.raw ) <- paste0( "A", seq_along( crx.raw ) )

mask <- apply( is.na( crx.raw ), 1, any )
crx <- crx.raw[ !mask, ]

## Normalize data
credit <- data.frame(
	  A1 = tolower( crx$A1 ) == 'a'
	, A2 = crx$A2
	, A3 = crx$A3
	, to_wide( crx$A4, "A4." )
	, to_wide( crx$A5, "A5." )
	, to_wide( crx$A6, "A6." )
	, to_wide( crx$A7, "A7." )
	, A8 = crx$A8
	, A9 = tolower( crx$A9 ) == 't'
	, A10 = tolower( crx$A10 ) == 't'
	, A11 = crx$A11
	, A12 = tolower( crx$A12 ) == 't'
	, to_wide( crx$A13, "A13." )
	, A14 = crx$A14
	, A15 = crx$A15
	, A16 = tolower( crx$A16 ) == '+' )

cbind( crx,  )

################################################################################
data(iris)
## Alias the data
d <- structure( iris, names = c( "sl", "sw", "pl", "pw", "tax" ) )

## Histogram of Sepal Length and Width
op <- par( mfcol = c( 1, 2 ) )
hist( d$sl, n = 20, xlab = "Sepal Length" ) ; hist( d$sw, n = 20, xlab = "Sepal Width" )
par( op )

## Define new bins for Sepal Length
a_sl <- c( 5.3, 6, 7.1 )

## Define new bins for Sepal Width
a_sw <- c( 2.5, 3.0, 3.6 )

## Split into bins
d <- cbind( d
	, sl.g = .bincode( d$sl, c( -Inf, a_sl, +Inf ) )
	, sw.g = .bincode( d$sw, c( -Inf, a_sw, +Inf ) ) )

## Construct a contingency table
ct_r <- tapply( d$sw, list( d$sl.g, d$sw.g ), length )
ct_l <- cbind( ct_r, rowSums( ct_r, na.rm = TRUE ) )
ct <- rbind( ct_l, colSums( ct_l, na.rm = TRUE ) )

## Display it!
image( ct )

################################################################################
data( iris )

x <- iris$Petal.Width
y <- iris$Petal.Length

bstrap <- function( y, x, H, N = 100 ) {
## Setup a data frmae for bootstrapping linear regression
	YX <- model.frame( y ~ x )
## Replicate the bootrstap sampling the given number of times
	replicate( n = N, {
## Draw a uniform sample with replacement
		bi <- sample( nrow( YX ), replace = TRUE )
## Estimate the model
		fi <- lm( y ~ x, data = YX[ bi, , drop = FALSE ] )
## Compute the bootstrapped t-statistic
		( coef( fi ) - H ) / sqrt( diag( vcov( fi ) ) )
	} )
}

bstrap( y, x, c(0,0), N=100)

################################################################################
n <- 4000
y <- 250

# set.seed(1)

ID <- sample( as.vector( sapply( LETTERS, paste0, LETTERS ) ), y )
d <- data.frame(
	p = sample( ID, n, replace = TRUE )
	, v = exp( rnorm( n ) * .1 )
	, t = sample( seq( as.Date( "2000-01-01" ), by = "day", len = n ), replace = TRUE ) )

	r0 <- do.call(rbind, lapply(split(d, d$p), function(d){
		d <- d[order(d$t),]
		d$v <- mapply(function(i,j) sum(d$v[i:j]),
			seq_along(d$t), .bincode(d$t+90, c(d$t,+Inf)))
		d[which.max(d$v),]
	}))

##############################################################################

library( mnormt )
r <- .6
dc <- function( u, v ) as.numeric(
	dmnorm(
		cbind( qnorm( u ), qnorm( v ) ),
		varcov = matrix( c( 1, r, r, 1 ), 2, 2 ) ) /
	dnorm( qnorm( u ) ) / dnorm( qnorm( v ) ) )

n <- 500
vectu <- seq( 1 / n, 1 - 1 / n, length = n - 1 )
matdc <- outer( vectu, vectu, dc )

contour( vectu, vectu, matdc, levels = c( .325, .944, 1.212, 1.250, 1.290, 1.656, 3.85 ), lwd = 2 )

require( ggplot2 )

DF <- data.frame(month = factor(month.abb, levels = month.abb),
                   freq = rpois(12, 80))
ggplot( DF, aes( x = month, y = freq ) ) +
	theme_bw( ) +
	coord_polar( ) +
	geom_bar( stat = "identity", fill = "pink2" )

################################################################################

# Attribution:
# Dubin, Robin A. (1992). Spatial autocorrelation and neighborhood quality. Regional Science and Urban Economics 22(3), 433-452.
ogrListLayers( "./baltim/baltim.shp" )
baltim <- readOGR(
	  "./baltim/baltim.shp"
	, "baltim", verbose = TRUE )

plot( baltim )

invisible( lapply( baltim@polygons, function( poly ) {
	polygon( poly@Polygons[[ 1 ]]@coords, col = "red", border = NA )
} ) )

xdbf<- read.dbf("./baltim/baltim.dbf")



ogrListLayers( "./QUAD250/quad250.shp" )

quad250 <- readOGR(
	  "./QUAD250/quad250.shp"
	, "quad250", verbose = TRUE )



library( maps )
library( rgdal )
library( foreign )
setwd( "/users/user/Desktop/studies 2014-2015/MIRKIN/datasets/" )
oz9799 <- readOGR(
	  "./oz9799/oz9799.shp"
	, "oz9799", verbose = TRUE )

plot(oz9799)


################################################################################
# Compuational intelligence
# Nonlinear regression using natrue-inspired apporach

X <- runif( 50 )
Y <- pmax( 2 * X^1.07 + rnorm( 50 ) * 2, 1.01 )


optim(c(0,1), function( p ) sum( (Y - p[1]*X^p[2])^2 ) )


## Determine the admissible area
step01 <- function( X, Y ) {
	ln_X <- log( X ) ; ln_Y <- log( Y )
	XY <- mapply(list, x=ln_X, y=ln_Y, SIMPLIFY = FALSE )
	# XY <- data.frame( x = ln_X, y = ln_Y )

	ZZ <- lapply( lapply( rev( seq_along( XY ) ), seq_len ), function( group ) {
		if( length( group ) < 2 ) return( NA )
		xy0 <- XY[ tail( group, n =  1 ) ][[1]]
		xy1 <- XY[ head( group, n = -1 ) ]
## Compute the slope
		coef <- sapply( xy1, function( pt ) {
			beta <- ( pt$y - xy0$y ) / ( pt$x - xy0$x )
			alpha <- pt$y -beta * pt$x
			c(alpha, beta)
		} )
		t(apply( coef, 1, range ))
	} )

	alpha <- sapply( ZZ, `[`, 1:2 ) ; beta <- sapply( ZZ, `[`, 3:4 )

	grid <- seq_along( XY )
	tapply( rep( XY, grid), rep(seq_along( XY ), rev( grid ) ), function( x ) {
		
		sapply( x[-1], function( y ) {
## Compute the slope
			
		} ) 
	} )


	b <- outer( ln_Y, ln_Y, `-` ) / outer( ln_X, ln_X, `-` )
	range( b, na.rm = TRUE )

}

step03 <- function( x, fn ) {
## Evaluate the 

}


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

## Do the PCA analysis
library( plyr )
## 1. Normalize the data: subtract centrality divide by scale
center <- sapply( wdbc[-1], mean )
stddev <- sqrt( diag( var( wdbc[-1] ) ) )

wdbc_s <- mapply( function( x, a, s ) ( x - a ) / s,
	wdbc[-1], center, stddev )

## 2. Compute the covariance matrix and its decomposition
ei <- eigen( crossprod( wdbc_s ) / nrow( wdbc_s ) )
pc_sdev <- sqrt( ei$values )
names( ei$vectors ) <- paste0( "PC", seq_along( pc_sdev ) )

wdbc_pc <- cbind( diagnosis = wdbc$diagnosis,
	as.data.frame( wdbc_s %*% ei$vectors ) )

ggplot( data = wdbc_pc ) + theme0( ) +
	geom_point( aes( V1, V2, colour = diagnosis ) )



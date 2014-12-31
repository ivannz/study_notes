rm( list = ls( all.names = TRUE) ) ; invisible( gc( ) )

setwd( "/users/user/Desktop/studies 2014-2015/Robust methods/assign02/tex" )

Sweave( file = "../assign02.Stex" )

##  Kullback-Leibler divergence between \mu_1 and \mu_2, parametrized by \alpha \in \clo{0,1}

ln_p_pq <- function( x, theta ) {
	d1 <- ( 1 - theta ) + 2 * theta * x
	d2 <- ( 1 - theta ) + theta / sqrt( - pi * log( x ) )
	log( d1 / d2 ) * d1
}
i_pq <- function( theta ) integrate( ln_p_pq, 0, 1, theta = theta, stop.on.error = FALSE, subdivisions = 200 )$value

ln_q_qp <- function( x, theta ) {
	d1 <- ( 1 - theta ) + 2 * theta * x
	d2 <- ( 1 - theta ) + theta / sqrt( - pi * log( x ) )
	log( d2 / d1 ) * d2
}
i_qp <- function( theta ) integrate( ln_q_qp, 0, 1, theta = theta, stop.on.error = FALSE, subdivisions = 200 )$value

x <- seq( 0, 1, by = 0.001 )
data <- data.frame( x = x, y_1 = sapply( x, i_qp ), y_2 = sapply( x, i_pq ) )

plot( x, sapply( x, i_qp ), type = "l", col = "black", lwd = 2 )
lines( x, sapply( x, i_pq ), type = "l", col = "red", lwd = 2 )



library( ggplot2 )
theme0 <-
  function( ... ) theme_bw( ) +
    theme_minimal( base_size = 18 ) # + theme( legend.position = "none" )

ggplot( data ) + theme0( ) +
  geom_line( aes( x = x, y = y_1, colour = "red" ) ) +
  geom_line( aes( x = x, y = y_2, colour = "blue" ) ) +
  labs( list( x = expression( theta ), y = NULL ) ) +
  xlim( 0, 1 )


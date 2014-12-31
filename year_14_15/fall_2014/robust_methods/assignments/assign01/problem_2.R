rm( list = ls( all.names = TRUE) ) ; invisible( gc( ) )

setwd( "/users/user/Desktop/studies 2014-2015/Robust methods/assign01/tex" )

Sweave( file = "../assign01.Stex" )

B <- pnorm( .5 )
C <- 8 * ( 1 - B ) / ( 3 * pi )

A <- c( C*pi/2, C*pi/2 + 2*Phi_5 - 1, C*pi/2 + 2*Phi_5 - 1 + C*pi/4 )

F <- function( x ) suppressWarnings(
  C * pi / 2 + ifelse( x < 0, C * atan( x ),
    B - 1 + ifelse( x <= 4, pnorm( ( x - 2 ) / 4 ),
      B + C/2 * atan( ( x - 4 )^2 ) ) ) )

F_inv <- function( u ) suppressWarnings(
  ifelse( u < C * pi / 2, tan( u / C - pi / 2 ),
    ifelse( u <= C * pi / 2 + 2 * B - 1,
      4 * qnorm( u - C * pi / 2 + 1 - B ) + 2,
        4 + sqrt( tan( 2 * ( u - C * pi / 2 - 2 * B + 1 ) / C ) ) ) ) )

u <- runif( 1000000 )
v <- F( e <- F_inv( u ) )
# qplot( e[e > -10 & e < 10], geom = "histogram", binwidth = 10 ) + theme_bw()

hist( e[e>-10 & e<10], n = 200 )


z <- matrix( rnorm( 3*100000 ), ncol=3 )
x <- c( 0, 0, 0 ) + chol( matrix( c( 16, 0, 0, 0, 9, 0, 0, 0, 9 ), 3 ) ) %*% z

x <- rmvnorm( n = 100000,
	mean = c( 0, 0, 0 ),
	sigma = matrix( c(
		   16, 0,  0
		,   0, 9,  0
		,   0, 0,  9 ), 3 ) )

xi <- ( 3 * x[ , 2 ] + x[ , 1 ] * x[ , 3 ] ) / sqrt( 9 + x[ , 1 ]^2 )

hist( xi / 3, n = 200, freq  =FALSE )
curve( dnorm, add = TRUE, n = 100, -5, 5, col = "red" )


## c\int_4^x \frac{s-4}{\brac{s-4}^4+1} ds 
## \frac{c}{2} \int_0^{x-4} \frac{2t}{t^4+1} dt 
## \frac{c}{2} \int_0^{{(x-4)}^2} \frac{1}{v^2+1} dv 
 

# \int_0^\infty \frac{t}{t^4 + 1} dt 
#  &= \obj{ x=t^2 } &=\frac{1}{2} \int_0^\infty \frac{1}{x^2+1} dx
#  &= \frac{1}{2} \frac{\pi}{2} = \frac{\pi}{4}

# \int_{-\infty}^{+\infty} f_X(x) dx 
#  = \int_0^4 \frac{1}{\sqrt{2\pi}\sigma} e^{-\frac{\brac{x-a}^2}{2\sigma^2}} dx
#   + c\int_{-\infty}^0 \frac{1}{x^2+1}dx  + c\int_4^{+\infty} \frac{x-4}{\brac{x-4}^4+1} dx
#  = \int_{\frac{-a}{\sigma}}^{\frac{4-a}{\sigma}} \frac{1}{\sqrt{2\pi}} e^{-\frac{t^2}{2}} dt +
#   + c\int_0^\infty \frac{1}{x^2+1}dx + c \int_0^\infty \frac{s}{s^4+1} ds
#  = \Phi\brac{\frac{4-a}{\sigma}}-\Phi\brac{\frac{-a}{\sigma}} +
#   + c\int_0^\infty \frac{1}{x^2+1}dx + \frac{c}{2} \int_0^\infty \frac{1}{v^2+1} dv
#  = \Phi\brac{\frac{4-a}{\sigma}}-\Phi\brac{\frac{-a}{\sigma}} +
#   + \frac{3 c}{2} \frac{\pi}{2} = 1

# c = \brac{\Phi(.5) - \Phi(-.5)} \frac{4}{3 \pi}
factors <-
  expand.grid( x1 = c(-1, 0, 1), x2 = c(1, 2) )[
    rep( 1 : 6, c( 2, 3, 1, 3, 1, 2 ) ), ]

X <- model.matrix( ~ 1 + x1 + x2 + x1*x2, factors )

C <- c( 0.3, -0.5, -0.5, 0.05 )

eta <- rnorm( nrow( X ) )

Y <- X %*% C + eta

## OLS
beta <- solve( crossprod( X ) ) %*% crossprod( X, Y )

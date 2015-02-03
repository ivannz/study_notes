import numpy as np
import numpy.random as rnd

import numpy.linalg as la

def synthfbmcircul( N, H, seed = None ) :
## Set both the values and the time scales
	sigma = 1.0 ; tmax = 1.0 ;
## The number of steps
	n = np.array( xrange( N ), dtype = np.dtype( "float64" ) )
## The base time increment of the mesh
	dt = tmax / N
## "Synthese de la covariance du fGn",
##  Synthesise the covariance of the fractional Brownian Motion
## Generate the first row of the 2Mx2M Toeplitz matrix
	H2 = 2 * H
	L = sigma / 2 * (
		np.abs( n - 1 ) ** H2 +
		np.abs( n + 1 ) ** H2 -
		2 * np.abs( n ) ** H2 ) * ( dt ** H2 ) ;
	print( "1" )
## Compute the FFT on real input (faster)
	L = np.real( np.fft.fft( np.append( L, L[1:-1][::-1] ) ) )
	print( "2" )
	w = rnd.randn( 2 * N - 2 ) * 1j + rnd.randn( 2 * N - 2 ) ;
# %% ATTENTION: ne pas utiliser ifft, qui utilise une normalisation differente
	print( "3" )
	z = np.fft.fft( np.sqrt( L / ( 2 * N - 2 ) ) * w )
## Prepend a zero to Z and get a cumulative sum
	print( "4" )
	z = np.append( 0, np.cumsum( z[ :N ] ) )
# % Processus increment normalise
## Dietrich, Newsam 1997: the real and imaginary parts of any m+1
##  consequtive entries yield two independent realizations of an FBM
	print( "5" )
	u = np.real( z[0:(N+1)] ) ; v = np.imag( z[0:(N+1)] )
	print( "6" )
	return( [np.real( z ), np.imag( z )] )
## Check the algorithm output and compare wit hthat of the original code

import numpy as np
import numpy.random as rnd
from numpy.fft import fft

def synthfbmcircul( N, H, seed = None ) :
	if seed is not None:
		rnd.seed( seed )
## Set both the values and the time scales
	sigma = 1.0 ; tmax = 1.0 ;
## The number of steps
	n = np.arange( N, dtype = np.float64 )
## The base time increment of the mesh
	dt = tmax / N
## "Synthese de la covariance du fGn",
##  Synthesise the covariance of the fractional Brownian Motion
## Generate the first row of the 2Mx2M Toeplitz matrix
	H2 = 2 * H
	L = sigma / 2 * ( dt ** H2 ) * (
		np.abs( n - 1 ) ** H2 + np.abs( n + 1 ) ** H2 - 2 * np.abs( n ) ** H2 ) ;
## Compute the FFT on real input (faster)
	L = np.real( fft( np.append( L, L[1:-1][::-1] ) ) )
	w = rnd.randn( 2 * N - 2 ) * 1j + rnd.randn( 2 * N - 2 ) ;
	# %% ATTENTION: ne pas utiliser ifft, qui utilise une normalisation differente
	z = fft( np.sqrt( L / ( 2 * N - 2 ) ) * w )
## Prepend a zero to Z and get a cumulative sum
	z = np.append( 0, np.cumsum( z[ :(N-1) ] ) )
## Dietrich, Newsam 1997: the real and imaginary parts of any m+1
##  consequtive entries yield two independent realizations of an FBM
# % Processus increment normalise
	return (np.real( z ), np.imag( z ))

# for mc in xrange( 10 ) :
# 	u, v = synthfbmcircul( 1024, .5 )
# 
# plt.plot(u, v, "b-", linewidth = 2 )
# plt.show()
# 
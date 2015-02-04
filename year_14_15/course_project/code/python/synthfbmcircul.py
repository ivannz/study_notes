## Fractional Brownian motion synthesis with circulant matrix method
##  of Dietrich, Newsam 1997.
## B(t) is synthetized for t in [0,tmax] and variance (of white gaussian noise) 
## equal to sigma2.
##
## Inputs:
##          N : number of samples
##          H : Hurst parameter (0 < H < 1)
##     [seed] : states (integers) for random number
##              generator (real sequence and imaginary sequence)
## Outputs:
##          B : H-fBm N-length trace
##          x : H-fGn : x[n] = B[n+1] - B[n]
##          w : white gaussian noise generator (complex)
##
##  >  > >>  Take N = 2^k + 1  << <  <
##
## Example :
##-----------
## N = 2^10+1 ;   % N=1025
## H = 1/3 ;
## [B,x,w] = synthfbmcircul(N,H) ;

## u = synthfbmcircul( 2**10 + 1, .3 )
## plt.plot(u[0], "b-", linewidth = 2 )
## plt.show()
## 

import numpy as np
import numpy.random as rnd
from numpy.fft import fft

def synthfbmcircul( N, H, sigma2 = 1.0, tmax = 1.0, seed = None, integrate = True ) :
	if seed is not None : rnd.seed( seed )
## The number of steps
	n = np.arange( N, dtype = np.float64 )
## The base time increment of the mesh
	dt = tmax / N
## "Synthese de la covariance du fGn",
##  Synthesise the covariance of the fractional Brownian Motion
## Generate the first row of the 2Mx2M Toeplitz matrix
	H2 = 2 * H
	L = sigma2 / 2 * ( dt ** H2 ) * (
		np.abs( n - 1 ) ** H2 + np.abs( n + 1 ) ** H2 - 2 * np.abs( n ) ** H2 ) ;
## Compute the FFT on real input (faster)
	L = np.real( fft( np.append( L, L[1:-1][::-1] ) ) )
	w = rnd.randn( 2 * N - 2 ) * 1j + rnd.randn( 2 * N - 2 ) ;
## Dietrich, Newsam 1997: the real and imaginary parts of any m+1
##  consequtive entries yield two independent realizations of an FBM
# %% ATTENTION: ne pas utiliser ifft, qui utilise une normalisation differente
	z = fft( np.sqrt( L / ( 2 * N - 2 ) ) * w )
## Prepend a zero to Z and get a cumulative sum
	if integrate :
# % Processus increment normalise
		z = np.append( 0, np.cumsum( z[ :N-1 ] ) )
	return ( np.real( z[ :N ] ), np.imag( z[ :N ] ) )

import matplotlib.pyplot as plt
u, v = synthfbmcircul( 2**10 + 1, .9, integrate = True, seed = 1 )
plt.plot(u, v, "r-", linewidth = 2 )
plt.show()

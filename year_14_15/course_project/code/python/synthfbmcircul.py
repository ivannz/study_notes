import numpy as np
import numpy.random as rnd
from numpy.fft import fft

## Genration of the fractional Gaussian noise with theh circulant matrix method
##  of Dietrich, Newsam 1997.
## Inputs:  N : number of samples
##          H : Hurst parameter (0 < H < 1)
##     [seed] : the initial state on the random number generator (for
## 				real and imaginary white noise sequences)

## def cplx_gaussian( M ) :
## 	for i in xrange( M ) :
## 		yield rnd.randn() + rnd.randn() * 1j
#### For large arrays use. However the fft fails: maybe use rfft?
##	 w = np.fromiter( cplx_gaussian( 2 * N - 2 ), dtype = np.complex128, count = 2 * N - 2 )

def synthgausscircul( N, H, variance = 1.0, seed = None ) :
## Begin with generation of the base white gaussian noise.
	if seed is not None : rnd.seed( seed )
	w = rnd.randn( 2 * N - 2 ) + rnd.randn( 2 * N - 2 ) * 1j
## "Synthese de la covariance du fGn", Synthesise the covariance
##  of the fractional Gaussian noise. This autocorrelation function
##  models long range (epochal) dependence.
## Generate the first row of the 2Mx2M Toeplitz matrix, where 2M = N + N-2
	n = np.arange( N, dtype = np.float64 )
	L = variance * ( np.abs( n - 1 ) ** (2 * H) + np.abs( n + 1 ) ** (2 * H) - 2 * np.abs( n ) ** (2 * H) ) / 2 ;
## Compute the convolution of the circulant row (of autocorrelations) with
##  some gaussian white noise
	L = np.real( fft( np.append( L, L[1:-1][::-1] ) ) )
# %% ATTENTION: ne pas utiliser ifft, qui utilise une normalisation differente
## F \times {(\tfrac{1}{2M}\Lambda)}^\tfrac{1}{2} \times w (see p.~1091)
	z = fft( np.sqrt( L / ( 2 * N - 2 ) ) * w )
## Dietrich, Newsam 1997: In our case the real and imaginary parts of any N
##  consequtive entries yield two independent realizations of \mathcal{N}_N(0,R)
##  where $R$ is the autocorrelation structure of an fBM
	return ( np.real( z[ :N ] ), np.imag( z[ :N ] ) )

## Fractional Brownian motion synthesis
## B(t) is synthetized for t in [0,tmax] and variance (of white gaussian noise) 
## equal to sigma2.
def synthfbmcircul( N, H, sigma2 = 1.0, tmax = 1.0, seed = None ) :
	dt = tmax / N
## for fBM use sigma2 = S dt^{2H}, with dt = (t_1-t_0)/N
	u, v = synthgausscircul( N, H, variance = sigma2 * ( dt * dt )**H, seed = seed )
## Return the integrated paths
# % Processus increment normalise
	return ( np.append( 0, np.cumsum( u[ :N-1 ] ) ), np.append( 0, np.cumsum( v[ :N-1 ] ) ) )

import matplotlib.pyplot as plt
u, v = synthgausscircul( 2**20 + 1, .85 )
u, v = synthfbmcircul( 2**20 + 1, .85 )
plt.plot(u, "r-", linewidth = 2 )
plt.plot(v, "b-", linewidth = 2 )
plt.show( )


## Example :
##  >  > >>  Take N = 2^k + 1  << <  <
##-----------
## N = 2^10+1 ;   % N=1025
## H = 1/3 ;
## [B,x,w] = synthfbmcircul(N,H) ;

## u = synthfbmcircul( 2**10 + 1, .3 )
## plt.plot(u[0], "b-", linewidth = 2 )
## plt.show()
## 

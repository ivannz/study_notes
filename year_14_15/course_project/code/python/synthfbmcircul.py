# -*- coding: UTF-8 -*-
import numpy as np

import numpy.random as rnd

from numpy.fft import fft, rfft

## Genration of the fractional Gaussian noise with theh circulant matrix method
##  of Dietrich, Newsam 1997.
## Inputs:  N : number of samples
##          H : Hurst parameter (0 < H < 1)
##     [seed] : the initial state on the random number generator (for
## 				real and imaginary white noise sequences)
class synth_fgn(object):
	"""A class to conveniently generate fractional Gaussian process using
		a circulant matrix method suggested by Dietrich and Newsam (1997)"""
	def __init__(self, N, H, sigma = 1.0 ) :
## Remeber the sample size
		self.__N = N
## The autocorrelation structure for the fBM is constant provided
##  the Hurst exponent and the size sample are fixed
## "Synthese de la covariance du fGn", Synthesise the covariance
##  of the fractional Gaussian noise. This autocorrelation function
##  models long range (epochal) dependence.
		R = np.arange( N, dtype = np.float64 )
## The noise autocorrelation structure is directly derivable from
##  the autocorrelation of the fBM:
##  r(s,t) = .5 * ( |s|^{2H}+|t|^{2H}-|s-t|^{2H} )
## If the noise is generated for an equally spaced. sampling of 
##  an fBM like process, then the autocorrelation function must
##  be multiplied by |Delta^{2H}
## Since FT is a linear transform (even the discrete one), this routine
##  can just generate a unit variance fractional gaussian noise.
		R = sigma * sigma * .5 * ( np.abs( R - 1 ) ** (2.0 * H) + np.abs( R + 1 ) ** (2.0 * H) - 2 * np.abs( R ) ** (2.0 * H) )
## Generate the first row of the 2Mx2M Toeplitz matrix, where 2M = N + N-2:
##  it should be r_0, ..., r_{N-1}, r_{N-2}, ..., r_1
		Z = np.real( rfft( np.append( R, R[::-1][1:-1] ) ) ) / ( 2 * N - 2 )
## The circulant matrix, defined by the autocorrelation structure above
##  is necessarily positive definite, which is equivalent to the FFT of any 
##  its row being nonegative. Due to numerical we truncate the close to zero
##  neagtive fourier coefficients.
		Z = np.sqrt( np.maximum( Z, 0.0, out = Z ) )
## Collect the frequencies of the specialized real FFT output by concatenating
##  F with F[::-1][1:-1] so that it mathces the output of the complex FFT. Real
##  FFT returns (2*N-2)/2+1 = N
		self.__acf_ft = np.append( Z, Z[::-1][1:-1] )
## The circulant embedding method actually generates a pair of independent
##  long-range dependent processes.
		self.__queue = []
## Use call semantics for brevity
	def __call__( self, seed = None ) :
## Generate the next sample only if needed.
		if not self.__queue :
			self.__queue.extend( self.__gen( seed ) )
## Return a pregenerated sample
		return self.__queue.pop( )
## The actual circulant generator
	def __gen( self, seed = None ) :
## Basically the idea is to utilize the convolution property of the Fourier
##  Transform and multiply the transform of the autocorrelation frunction by
##  the independent gaussian white noise in the frequency domain and then
##  get back to the time domain.
## cf. http://www.thefouriertransform.com/transform/properties.php
		if seed is not None : rnd.seed( seed )
## Begin with generation of the standarized white gaussian noise.
		W = rnd.randn( 2 * self.__N - 2 ) + rnd.randn( 2 * self.__N - 2 ) * 1j
## Compute the convolution of the circulant row (of autocorrelations) with
##  some gaussian white noise.
# %% ATTENTION: ne pas utiliser ifft, qui utilise une normalisation differente
## F \times {(\tfrac{1}{2M}\Lambda)}^\tfrac{1}{2} \times w (see p.~1091 [Dietrich,Newsam; 1997])
		W = fft( self.__acf_ft * W )
## Dietrich, Newsam 1997: "In our case the real and imaginary parts of any N
##  consequtive entries yield two independent realizations of \mathcal{N}_N(0,R)
##  where $R$ is the autocorrelation structure of an fBM."
## Therefore take the first N complex draws.
		return ( np.real( W[ :self.__N ] ), np.imag( W[ :self.__N ] ) )
	def reset( self ):
		del self.__queue[:]

## A class for generating fractional brownian motion
## Actually fBM is a Hermite process of order 1
class synth_fbm( synth_fgn ):
	"""Generator: Produce sample paths of a fractional Brownian Motion with
	the specified Hurst exponent"""
## For better performance N should be a power of two.
	def __init__(self, N, H, tmax = 1.0 ):
		dt = tmax / N
		self.__t = tmax * np.arange( N, dtype = np.float ) / ( N - 1 )
		super(synth_fbm, self).__init__(N + 1, H, dt ** H)
	def __call__( self, seed = None ) :
		increments = super(synth_fbm, self).__call__( seed )
		return self.__t, np.cumsum( increments[:-1] )
	def reset( self ):
		super(synth_fbm, self).reset( )

# def synthgausscircul( N, H, variance = 1.0, seed = None ) :
# 	if seed is not None : rnd.seed( seed )
# 	w = rnd.randn( 2 * N - 2 ) + rnd.randn( 2 * N - 2 ) * 1j
# 	n = np.arange( N, dtype = np.float64 )
# 	L = variance * ( np.abs( n - 1 ) ** (2.0 * H) + np.abs( n + 1 ) ** (2.0 * H) - 2 * np.abs( n ) ** (2.0 * H) ) / 2 ;
# 	Z = np.sqrt( np.real( rfft( np.append( L, L[::-1][1:-1] ) ) ) / ( 2 * N - 2 ) )
# 	Z = np.append( Z, Z[::-1][1:-1] )
# 	Z = fft( Z * w )
# 	return ( np.real( Z[ :N ] ), np.imag( Z[ :N ] ) )
# def synthfbmcircul( N, H, sigma2 = 1.0, tmax = 1.0, seed = None ) :
# 	dt = tmax / N
# 	u, v = synthgausscircul( N, H, variance = sigma2 * ( dt * dt ) ** H, seed = seed )
# 	return ( np.append( 0, np.cumsum( u[ :N-1 ] ) ), np.append( 0, np.cumsum( v[ :N-1 ] ) ) )

# Example
##  >  > >>  Take N = 2^k + 1  << <  <for better performance
# import matplotlib.pyplot as plt
# gen = synth_fgn( 2**10 + 1, .995 )
# u, v = gen(), gen()
# plt.plot(u, "r-", linewidth = 2 )
# plt.plot(v, "b-", linewidth = 2 )
# plt.show( )


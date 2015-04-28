# -*- coding: UTF-8 -*-
import numpy as np
from numpy.fft import rfft, fft

class fgn( object ) :
	"""A class to generate fractional Gaussian process of fixed length using
	a circulant matrix method suggested by Dietrich and Newsam (1997). For the
	best performance N-1 should be a power of two."""
## For better performance N-1 should be a power of two.
	def __init__( self, N, H = .5, sigma = 1.0 ) :
## The autocorrelation structure for the fBM is constant provided the Hurst exponent
##  and the size sample are fixed. "Synthese de la covariance du fGn", Synthesise
##  the covariance of the fractional Gaussian noise. This autocorrelation function
##  models long range (epochal) dependence.
		R = np.arange( N, dtype = np.float64 )
## The noise autocorrelation structure is directly derivable from the autocorrelation
##  of the time-continuous fBM:
##     r(s,t) = .5 * ( |s|^{2H}+|t|^{2H}-|s-t|^{2H} )
## If the noise is generated for an equally spaced. sampling of an fBM like process,
##  then the autocorrelation function must be multiplied by ∆^{2H}. Since Fourier
##  Transform is linear (even the discrete one), this routine can just generate a unit
##  variance fractional Gaussian noise.
		R = sigma * sigma * .5 * (
			  np.abs( R - 1 ) ** ( 2.0 * H )
			+ np.abs( R + 1 ) ** ( 2.0 * H )
			- 2 * np.abs( R ) ** ( 2.0 * H ) )
## Generate the first row of the 2Mx2M Toeplitz matrix, where 2M = N + N-2: it should
##  be [ r_0, ..., r_{N-1}, r_{N-2}, ..., r_1 ]
		Z = np.real( rfft( np.append( R, R[::-1][1:-1] ) ) )
		del R
## The circulant matrix, defined by the autocorrelation structure above is necessarily
##  positive definite, which is equivalent to the FFT of any its row being non-negative.
## Due to numerical round-off errors we truncate close to zero negative Fourier
##  coefficients.
		Z = np.sqrt( np.maximum( Z, 0.0 ) / ( 2 * N - 2 ) )
## Collect the frequencies of the specialized real FFT output by concatenating Z with
##  Z[::-1][1:-1] so that it matches the output of the complex FFT. Real FFT returns
##  a vector of (2*N-2)/2+1 = N coefficients.
		self.__acf_ft = np.append( Z, Z[::-1][1:-1] )
		del Z
## The circulant embedding method actually generates a pair of independent long-range
##  dependent processes.
		self.__queue = list( )
## Remember the sample size
		self.__N = N
## Setup a local rng
		self.__np_rand = None
## fGn generator via circulant embedding method
	def __gen( self ) :
## Basically the idea is to utilize the convolution property of the Fourier Transform
##  and multiply the transform of the autocorrelation function by the independent
##  Gaussian white noise in the frequency domain and then get back to the time domain.
##    cf. \url{ http://www.thefouriertransform.com/transform/properties.php }
## Begin with generation of the Gaussian white noise with unit variance and zero mean.
		W = self.__np_rand.randn( 2 * self.__N - 2 ) + self.__np_rand.randn( 2 * self.__N - 2 ) * 1j
## Compute the convolution of the circulant row (of autocorrelations) with the noise.
## "%% ATTENTION: ne pas utiliser ifft, qui utilise une normalisation differente"
## Compute this (see p.~1091 [Dietrich, Newsam; 1997]) :
##  F \times (\frac{1}{2M}\Lambda)^\frac{1}{2} \times w
		W = fft( self.__acf_ft * W )
## [Dietrich, Newsam; 1997] write : "In our case the real and imaginary parts of any N
##  consecutive entries yield two independent realizations of \mathcal{N}_N(0,R) where
##  $R$ is the autocorrelation structure of an fBM."
##  Therefore take the first N complex draws to get a pair of independent realizations.
		return ( np.real( W[ :self.__N ] ), np.imag( W[ :self.__N ] ) )
## Reset the internal state of the generator
	def reset( self ) :
		del self.__queue[:]
## A visbile function, to generate the sample
	def __call__( self ) :
## Generate the next sample only if needed.
		if not self.__queue :
			self.__queue.extend( self.__gen( ) )
## Return a pregenerated sample
		return self.__queue.pop( )
	def set_rnd( self, numpy_random ) :
		self.reset( )
		self.__np_rand = numpy_random


class fbm( fgn ):
	"""A derived class to produce sample paths of a Fractional Brownian Motion with
	a specified fractional integration parameter (the Hurst exponent). For the best
	performance N-1 should be a power of two."""
	def __init__(self, N, H = 0.5 ) :
		self.__t = np.arange( N, dtype = np.float ) / ( N - 1 )
		fgn.__init__( self, N, H, sigma = ( 1.0 / N ) ** H )
	def __call__( self ) :
		increments = super( fbm, self ).__call__( )
		return self.__t, np.concatenate( ( [ 0 ], np.cumsum( increments[ : -1 ] ) ) )
	def reset( self ):
		super( fbm, self ).reset( )
	def set_rnd( self, numpy_random ) :
		super( fbm, self ).set_rnd( numpy_random )
		self.reset( )



class test_gen(object):
	def __init__( self, N ) :
		self.__queue = list( )
		self.__N = N
		self.__np_rand = None
	def set_rnd( self, numpy_random ) :
		self.__np_rand = numpy_random
	def __call__( self ) :
		if not self.__queue :
			self.__queue.extend( self.__gen( )[::-1] )
		return self.__queue.pop( )
	def reset( self ) :
		del self.__queue[ : ]
	def __gen( self ) :
		return [( 0, self.__np_rand.randn( 1 ), ) for i in range( self.__N ) ]
	def info( self ) :
		return ( id( self.__queue ), )

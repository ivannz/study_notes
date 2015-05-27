# -*- coding: UTF-8 -*-
import numpy as np

from fgn import fgn

from numpy.polynomial.hermite_e import hermeval

class hermite( fgn ) :
	"""A derived class to produce sample paths of a Hermite process of order d
	with a specified fractional integration parameter (the Hurst exponent). For
	the best performance N-1 should be a power of two."""
	def __init__(self, N, d = 2, H = 0.5, K = 16, time = False, **kwargs ) :
		self.__t = np.empty( 0, np.float ) if not time else np.arange( N, dtype = np.float ) / ( N - 1 )
## It is imperative that the varince of the fGn be unity
## The hurst exponent for this process is H = 1 + d * ( Hfgn - 1 )
		fgn.__init__( self, K * ( N - 1 ) + 1, H = ( H + d - 1.0 ) / d, sigma = 1.0, **kwargs )
## The downsampling parameter (K) is actually the index of the process,
##  which when it tends to infinity, converges in distribution to the
##  Rosenblatt process (or in general to a Hermite process).
## The non-central limit theorem:
##  Z^k(t) = \frac{1}{n^\alpha}\sum_{j=1}^{\left\lfloor kt\right\rfloor} H(\xi_j)
## Converges to $Z_\frac{\alpha}{2}(t)$ -- a hermite process 
## Increasing downsample gives better approximation. In theory it should
##  tend to infinity. This is a serious drawback.
## c.f. [Abry, Pipiras; 2005]
		self.__K = K
		self.__H = H
## Define the order of the Hermite polynomial
		self.__coef = np.zeros( d + 1, np.float )
		self.__coef[ d ] = 1
	def __call__( self ) :
## Generate values of a hermite polynomial of the given order at the values
##  of a fractional Gaussian Noise with the specified hurst index.
		increments = hermeval( fgn.__call__( self ), self.__coef )
## The renorm-group transformation, without the renormalisation by the n^{-H}
		increments = np.cumsum( increments )[ self.__K-1::self.__K ] / ( self.__K ** self.__H )
		# return self.__t, np.concatenate( ( [ 0 ], increments ) )
		return self.__t, np.concatenate( ( [ 0 ], increments / np.max( increments ) ) )
	def reset( self ):
		super( hermite, self ).reset( )
	def set_rnd( self, numpy_random ) :
		super( hermite, self ).set_rnd( numpy_random )
		self.reset( )

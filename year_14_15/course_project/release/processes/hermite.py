# -*- coding: UTF-8 -*-
import numpy as np
from numpy.polynomial.hermite_e import hermeval

from processes.fgn_pyfftw import fgn
# from processes.fgn_numpy import fgn

# def hermite( fgn, N, D = 2, H = 0.5, K = 16, time = False, **kwargs ) :
	# return hermite( N, D, H, K, time, **kwargs )

class hermite( fgn ) :
	"""A derived class to produce sample paths of a Hermite process of order d
	with a specified fractional integration parameter (the Hurst exponent). For
	the best performance N-1 should be a power of two."""
	def __init__( self, N, D = 2, H = 0.5, K = 16, time = False, **kwargs ) :
		self.__t = np.empty( 0, np.float ) if not time else np.arange( N, dtype = np.float ) / ( N - 1 )
## It is imperative that the varince of the fGn be unity
## The hurst exponent for this process is H = 1 + D * ( Hfgn - 1 )
		fgn.__init__( self, K * ( N - 1 ) + 1, H = ( H + D - 1.0 ) / D, sigma = 1.0, **kwargs )
## The downsampling parameter (K) is actually the index of the process,
##  which when it tends to infinity, converges in distribution to the
##  Rosenblatt process (or in general to a Hermite process).
## The non-central limit theorem:
##  Z^k(t) = \frac{1}{n^\alpha}\sum_{j=1}^{\left\lfloor kt\right\rfloor} H(\xi_j)
## Converges to $Z_\frac{\alpha}{2}(t)$ -- a hermite process 
## Increasing downsample gives better approximation. In theory it should
##  tend to infinity. This is a serious drawback.
## c.f. [Abry, Pipiras; 2005]
		self.__K, self.__H = K, H
## Define the order of the Hermite polynomial
		self.__coef = np.zeros( D + 1, np.float )
		self.__coef[ D ] = 1
	def reset( self ):
		super( hermite, self ).reset( )
	def initialize( self, numpy_random_state ) :
		super( hermite, self ).initialize( numpy_random_state )
	def __call__( self ) :
## Generate values of a hermite polynomial of the given order at the values
##  of a fractional Gaussian Noise with the specified hurst index.
		increments = hermeval( super( hermite, self ).__call__( ), self.__coef )
## The renorm-group transformation, without the renormalisation by the n^{-H}
		increments = np.cumsum( increments )[ self.__K-1::self.__K ] ## / ( self.__K ** self.__H )
		return self.__t, np.concatenate( ( [ 0 ], increments / np.max( np.abs( increments ) ) ) )

#! /usr/bin/env python
# -*- coding: UTF-8 -*-
import numpy as np

from synthfbmcircul import synth_fgn as synth_fgn

## The downsample parameter (K) is actually the index of the process,
##  which when it tends to infinity, converges in distribution to
##  the Rosenblatt process (or in general to a Hermite process).
## The non-central limit theorem:
##  Z^k(t) = \frac{1}{n^\alpha}\sum_{j=1}^{\left\lfloor kt\right\rfloor} H(\xi_j)
## Converges to $Z_\frac{\alpha}{2}(t)$ -- a hermite process 
## Increasing downsample gives better approximation. In theory
##  it should tend to infinity. This is a serious drawback.
## C.f. [Abry, Pipiras; 2005]

## Hermite process of order 2
class synth_Rosenblatt( synth_fgn ):
	def __init__(self, N, H, K = 16 ):
		self.__t = np.arange( N, dtype = np.float ) / ( N - 1 )
## It is imperative that the vraince of the fGn be unity
		super(synth_Rosenblatt, self).__init__( N * K + 1, ( H + 1.0 ) / 2.0, 1.0 )
## The hurst exponent for this process is H = 1 + 2 * ( Hfgn - 1 )
		self.__K = K
	def __call__( self, seed = None ) :
## Generate fractional Gaussian Noise
		u = super(synth_Rosenblatt, self).__call__( seed )
## The renorm-group transformation, without the renormalisation by the n^{-H}
		u = np.cumsum( u*u - 1 )[ (self.__K-1): :self.__K ]
		return self.__t, u / np.max( np.abs( u ) )
	def reset( self ):
		super(synth_Rosenblatt, self).reset( )

## Hermite process of order 3
class synth_Hermite3( synth_fgn ):
	def __init__(self, N, H, K = 16 ):
		self.__t = np.arange( N, dtype = np.float ) / ( N - 1 )
## It is imperative that the vraince of the fGn be unity
		super(synth_Hermite3, self).__init__( N * K + 1, ( H + 2.0 ) / 3.0, 1.0 )
## The hurst exponent for this process is H = 1 + 3 * ( Hfgn - 1 )
		self.__K = K
	def __call__( self, seed = None ) :
## Generate fractional Gaussian Noise
		u = super(synth_Hermite3, self).__call__( seed )
## The renorm-group transformation, without the renormalisation by the n^{-H}
		u = np.cumsum( ( u*u - 3 ) * u )[ (self.__K-1): :self.__K ]
		return self.__t, u / np.max( np.abs( u ) )
	def reset( self ):
		super(synth_Hermite3, self).reset( )

## Hermite process of order 3
class synth_Hermite4( synth_fgn ):
	def __init__(self, N, H, K = 16 ):
		self.__t = np.arange( N, dtype = np.float ) / ( N - 1 )
## It is imperative that the vraince of the fGn be unity
		super(synth_Hermite4, self).__init__( N * K + 1, ( H + 3.0 ) / 4.0, 1.0 )
## The hurst exponent for this process is H = 1 + 4 * ( Hfgn - 1 )
		self.__K = K
	def __call__( self, seed = None ) :
## Generate a stationary Gaussian process with long-range dependence:
##  the covariance \mathbb{E}(\xi_0\xi_k) = k^{2H} L(k)
##  where L(k) is a \textbf{s}lowly \textbf{V}arying \textbf{F}unction,
##  in our case a constant. An SVF has \lim_{s\to\infty} \frac{L(st)}{L(s)} = const 
		u = super(synth_Hermite4, self).__call__( seed )
## Transform the generated process by a Hermite polynomian of order 4
## Square the generated values
		u = u*u - 3.0
## The synthesized process should have unit variance at t=1. This
##  normalisation does not guarantee this.
## The renorm-group transformation, without the renormalisation by the n^{-H}
		u = np.cumsum( u*u - 6.0 )[ (self.__K-1): :self.__K ]
		return self.__t, u / np.max( np.abs( u ) )
	def reset( self ):
		super(synth_Hermite4, self).reset( )

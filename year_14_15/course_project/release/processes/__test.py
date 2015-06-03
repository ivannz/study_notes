# -*- coding: UTF-8 -*-
import numpy as np

class test_gen( object ) :
	def __init__( self, N, **kwargs ) :
		self.__cache = list( )
		self.__N = N
		self.__np_rand = None
	def initialize( self ) :
		pass
	def set_rnd( self, numpy_random ) :
		self.__np_rand = numpy_random
	def __call__( self ) :
		if not self.__cache :
			self.__cache.extend( self.__gen( )[::-1] )
		return self.__cache.pop( )
	def reset( self ) :
		del self.__cache[ : ]
	def __gen( self ) :
		return [( 0, self.__np_rand.randn( 1 ), ) for i in range( self.__N ) ]
	def info( self ) :
		return ( id( self.__cache ), )

from numpy.random import RandomState

from processes.fgn_numpy import fgn as fgn_numpy
from processes.fgn_pyfftw import fgn as fgn_fftw

from processes.fbm import fbm
from processes.hermite import hermite

## Test FGN
N, H, sigma = 2**21 + 1, 0.76, 1.0

g1 = fgn_numpy( N, H, sigma = 1.0 )
g2 = fgn_fftw( N, H, sigma = 1.0 )

g1.initialize( RandomState( 321 ) )
g2.initialize( RandomState( 321 ) )

x1, x2 = g1( ), g2( )
assert( np.allclose( x1, x2 ) )
del x1, x2
del g1, g2

## Test hermite
N, H, D, K = 2**16 + 1, 0.76, 2, 16
g1 = hermite( fgn_numpy, N, D, H, K, time = True )
g2 = hermite( fgn_fftw, N, D, H, K, time = True )

g1.initialize( RandomState( 123 ) )
g2.initialize( RandomState( 123 ) )

(t1, x1), (t2, x2) = g1( ), g2( )
assert( np.allclose( x1, x2 ) )
del t1, t2, x1, x2
del g1, g2

## Test FBM
N, H = 2**21 + 1, 0.76
g1 = fbm( fgn_numpy, N, H, time = True )
g2 = fbm( fgn_fftw, N, H, time = True )

g1.initialize( RandomState( 654 ) )
g2.initialize( RandomState( 654 ) )

(t1, x1), (t2, x2) = g1( ), g2( )
assert( np.allclose( x1, x2 ) )
del t1, t2, x1, x2
del g1, g2


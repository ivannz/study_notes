# -*- coding: UTF-8 -*-
import numpy as np

class empty( object ) :
	def __init__( self, *wargs, **kwargs ) :
		pass
	def set_rnd( self, *wargs, **kwargs ) :
		pass
	def __call__( self, *wargs, **kwargs ) :
		return None
	def reset( self, *wargs, **kwargs ) :
		pass

class test_gen( object ) :
	def __init__( self, N, **kwargs ) :
		self.__cache = list( )
		self.__N = N
		self.__np_rand = None
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


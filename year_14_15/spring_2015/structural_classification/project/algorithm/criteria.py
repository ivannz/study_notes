#! /usr/bin/env python
# -*- coding: UTF-8 -*-

import numpy as np
from scipy.spatial.distance import pdist, squareform, cdist
from numpy import argpartition as argtop

class criterion(object):
	"""Класс реализует оценку качества разбиения исходных данных.
	Параметры инициализации:
	-- $\\alpha$, $p$ параметры трансформации расстояния в меру близости:
		$$K(x,y) = e^{ -\\alpha d(x,y)^p }$$
	Использование:
		crit = partition_criterion(
			data = данные,
			alpha = масштаб,
			p = степень )
	## Критерий внутриклассовой близости
		crit.sim_within( pi = разбиение на классы в формате списка списков )
	## Критерий межклассовой близости
		crit.sim_between( pi = разбиение на классы в формате списка списков )
	## Комбинированный критерий тербующий похожести внутри класса и
	##  различий между классами
		crit( pi = разбиение на классы в формате списка списков )
	"""
	def __init__(self, data, fun, alpha = .5, p = 2, **kwargs):
		super(criterion, self).__init__( )
		print alpha, p
		sim_matrix = squareform( np.exp( - alpha * ( pdist( data, **kwargs ) ** p ) ) )
## Reomve the diagonal, since self-similarity should not contribute to
##  the clustering quality criteria.
		np.fill_diagonal( sim_matrix, 0 )
		self.__similarity = sim_matrix
		self.__data = data
		self.__fun = fun
	def evaluate( self, partition ) :
		return self.__fun( partition, self.__similarity, self.__data )
	def find_candidates( self, src, dst, s = 1 ) :
		if s >= len( src ) :
			return set( src )
## Compute the current centroid of the dst-class : O( n )
		centroid = np.mean( self.__data[ dst ], axis = 0 )
## Compute the distances to the centroid : O( n )
		distances = cdist( [ centroid ], self.__data[ src ] )
## Find s nearest neighbours using introselect ( min-heap + quickselect ): <= O( n log s) + O( s log s)
		return set( argtop( distances[ 0 ], s )[ :s ] )

def sim_ratio( partition, sim, *wargs, **kwargs ) :
	return ( sim_within( partition, sim, *wargs, **kwargs ) / 
		sim_between( partition, sim, *wargs, **kwargs ) )

def sim_diff( partition, sim, *wargs, **kwargs ) :
	return ( sim_within( partition, sim, *wargs, **kwargs ) - 10 * 
		sim_between( partition, sim, *wargs, **kwargs ) )

def sim_diff_ratio( partition, sim, *wargs, **kwargs ) :
	I1 = sim_within( partition, sim, *wargs, **kwargs )
	I2 = sim_between( partition, sim, *wargs, **kwargs )
	return ( I1 - 10 * I2 ) / ( I1 + 10 * I2 )

def sim_within( partition, sim, *wargs, **kwargs ) :
	# if np.any( [ len( i ) < 2 for i in partition ] ) : return 0 ##-np.inf
	return np.mean( [ __sim_AA( sim, I ) for I in partition ], dtype = np.float64 )

def sim_between( partition, sim, *wargs, **kwargs ) :
	from itertools import combinations as comb
	return np.mean( [ __sim_AB( sim, I, J ) for I, J in comb( partition, 2 ) ], dtype = np.float64 )

def __sim_AA( sim, A ) :
## Atomic or empty classes must be heavily penalized, but not so much as to 
##  inhibit any data point transfers.
	if len( A ) < 1 : return -np.inf
	if len( A ) < 2 : return 0
	within = np.mean( sim[ np.ix_( A, A ) ], dtype = np.float64 )
## Reomve the diagonal, since self-similarity should not contribute to
##  the clustering quality criteria.
	diag = np.mean( sim[ A, A ], dtype = np.float64 )
	return ( len( A ) * within - diag ) / ( len( A ) - 1.0 )

def __sim_AB( sim, A, B ) :
## Implement average linkage
	return np.mean( sim[ np.ix_( A, B ) ], dtype = np.float64 )

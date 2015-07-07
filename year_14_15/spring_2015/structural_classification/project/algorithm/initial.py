#! /usr/bin/env python
# -*- coding: UTF-8 -*-

import numpy as __np
from scipy.cluster.hierarchy import linkage as __linkage	

def empty( data, num_classes = 5, **kwargs ) :
	"""Пустое разбиение"""
	return [ list( ) for k in range( num_classes ) ]

def hierarchical( data, num_classes = 5, **kwargs ) :
	"""Функция для порождения начального разбиения методом агломеративной
	иерархической кластеризации.
	Суть метода: на исходных данных строится иерархия классов согласно мере близости,
	указанной при вызове. Построение производится снизу--вверх последовательным
	объедининем пар множеств разбиения. По завершении производится откат на 'num_classes'
	шагов назад до той итерации, на которой количество классов уменьшилось до
	'num_classes' штук. Существенный недостаток методоа состоит в том, что при выборе
	оптимистичного метода оценки схожести множеств (single linkage), возможно порождение
	атомарных классов, неподходящих для m-локальной оптимизации.
	Использование:
		pi0 = hierarchical(
			data = данные,
			num_classes = количество классов,
			параметры вызова для scipy.cluster.hierarchy.linkage(...) )
	"""
	clust = __linkage( data, **kwargs )
## Scalp off the last iterations of agglomerative clustering
##  to figure out the partition with num_classes classes
	C = __np.sort( clust[ -num_classes+1:, :2 ].ravel( ) )[ :num_classes ]
	return [ __hierarchical_rebuild( clust, c ) for c in C ]

## Rebuild the class synthesized at k-th iteration of
##   agglomerative hierarchical clustering
def __hierarchical_rebuild( clust, iter ) :
	if iter <= clust.shape[ 0 ] : return [ int( iter ) ]
	i, j = clust[ iter - clust.shape[ 0 ] - 1, :2 ]
	return __hierarchical_rebuild( clust, i ) + __hierarchical_rebuild( clust, j )

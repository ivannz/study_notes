# -*- coding: UTF-8 -*-
## Base modules

import scipy.sparse as sp
import numpy as np
from collections import deque

import pandas as pd

## Read a sparese adjacency matrix from a two-column CSV file
def __csr_from_csv( file_name, **kwargs ) :
	return __csr_from_pandas( pd.read_csv( file_name, **kwargs ) )

## Creates a sparse matrix from a two-column source-destination dataframe
def __csr_from_pandas( df, symmetrize = False ) :
	return __csr_from_endpoints( df.values[ :, 0 ],
		df.values[ :, 1 ], symmetrize = symmetrize )

def __csr_from_endpoints( u, v, symmetrize = False ) :
	assert( len( u ) == len( v ) )
## Convert to a COO matrix
	if not symmetrize :
		adj = sp.coo_matrix( ( np.ones( len( u ), dtype = np.float ), ( u, v ) ) )
	else :
		adj = sp.coo_matrix( ( np.ones( len( u ) + len( v ), dtype = np.float ),
			( np.concatenate( ( u, v ) ), np.concatenate( ( v, u ) )) ) )
## Convert to CSR and remove duplicates
	adj = adj.tocsr( ) ; adj.data[ : ] = 1
	return adj

def __sparse_bfs( A, sources, num_nodes = np.inf, max_hops = np.inf ) :
	sources = np.asarray( sources, np.int )
## Initialize the hops array
	dist = np.full( A.shape[ 0 ], np.inf, np.float )
## THe source is immediately reachable
	dist[ sources ] = 0.0
## Setup the vertex traversal schedule.
	Q = deque( sources )
## Setup the counter of visited nodes
	num_visited = 0
## If the allotted number of nodes has been exceeded, break the cycle.
	while Q :
## Get the current vertex from the top of the FIFO queue
		v = Q.popleft( )
##  ... find its nerighbours (A is CSR)
		N = A[ v, : ].nonzero( )[ 1 ]
##  ... keep those that were not visited
		N = N[ np.isinf( dist[ N ] ) ]
## Add the mto the queue
		if len( N ) > 0 :
			dist[ N ] = 1.0 + dist[ v ]
## Nodes farther than max_hops away from the sources are not traversed.
			if 1.0 + dist[ v ] < max_hops :
				Q.extend( N )
## Unless the current vertex is the source, increase the number of visited nodes.
		if dist[ v ] > 0 :
			num_visited += len( N )
		if num_visited >= num_nodes : 
			break
	return dist


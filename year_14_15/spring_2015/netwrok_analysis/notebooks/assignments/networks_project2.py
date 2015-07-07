import pandas as pd, numpy as np, scipy.sparse as sp
import os, gc, regex as re, time as tm

import matplotlib.pyplot as plt
# %matplotlib inline

DATADIR = os.path.realpath( os.path.join( ".", "data", "proj02" ) )

raw_dblp_file = os.path.join( DATADIR, "dblp_2000.csv.gz" )
cached_dblp_file = os.path.join( DATADIR, "dblp_2000.ppdf" )
cached_author_index = os.path.join( DATADIR, "dblp_2000_authors.txt" )

## Return a mask of elements of b found in a: optimal for numeric arrays
def match( a, v, return_indices = False ) :
	index = np.argsort( a )
## Get insertion indices
	sorted_index = np.searchsorted( a, v, sorter = index )
## Truncate the indices by the length of a
	index = np.take( index, sorted_index, mode = "clip" )
	mask = a[ index ] == v
## return
	if return_indices :
		return mask, index[ mask ]
	return mask

## Convert the edgelist into sparse matrix
def to_sparse_coo( u, v, shape, dtype = np.int32 ) :
## Create a COOrdinate sparse matrix from the given ij-indices
	assert( len( u ) == len( v ) )
	return sp.coo_matrix( (
			np.ones( len( u ) + len( v ), dtype = dtype ), (
				np.concatenate( ( u, v ) ), np.concatenate( ( v, u ) ) )
		), shape = shape )
## Remeber: when converting COO to CSR/CSC the duplicate coordinate
##  entries are summed!

## Create cache if necessary
tick = tm.time( )
if not os.path.exists( cached_dblp_file ) :
	## Load the csv file into a dataframe
	dblp = pd.read_csv( raw_dblp_file, # nrows = 10000,
	## On-the-fly decompression
			compression = "gzip", header = None, quoting = 0,
	## Assign column headers
			names = [ 'author1', 'author2', 'year', ], encoding = "utf-8" )
	## Finish
	tock = tm.time( )
	print "Raw DBLP read in %.3f sec." % ( tock - tick, )
## Map author names to ids
	from sklearn.preprocessing import LabelEncoder
	le = LabelEncoder( ).fit( np.concatenate( (
		dblp["author1"].values, dblp["author2"].values, ) ) )
	dblp_author_index = le.classes_
	for col in [ 'author1', 'author2', ] :
		dblp[ col ] = le.transform( dblp[ col ] )
## Cache
	dblp.to_pickle( cached_dblp_file )
	with open( cached_author_index, "w" ) as out :
		for label in le.classes_ :
			out.write( label.strip( ).encode( "utf-8" ) + "\n" )
	del dblp, le
## Finish
	tick = tm.time( )
	print "Preprocessing took %.3f sec." % ( tick - tock, )
else :
## Load the database from pickled format
	dblp = pd.read_pickle( cached_dblp_file )
## Read the dictionary of authors
	with open( cached_author_index, "r" ) as author_index :
		dblp_author_index = [ line.decode( "utf-8" ) for line in author_index ]

## Report
tock = tm.time( )
print "DBLP loaded in %.3f sec." % ( tock - tick, )


## 1. Prepare the train sample
## 1.1 First preprocess the pre 2010 data.
pre2010 = dblp[ dblp.year <= 2010 ].copy( )

## 1.2.1 Reencode the vertices in a less wasteful manner
from sklearn.preprocessing import LabelEncoder
le = LabelEncoder( ).fit( np.concatenate( ( pre2010[ "author1" ].values, pre2010[ "author2" ].values ) ) ) 
pre2010_values = le.classes_

## 1.2.2 Recode the edge data
for col in [ 'author1', 'author2', ] :
    pre2010[ col ] = le.transform( pre2010[ col ] )      

## 1.3.1 Create a symmetric adjacency matrix
pre2010_adj = to_sparse_coo(
    pre2010[ "author1" ].values, pre2010[ "author2" ].values,
    shape = 2 * [ len( le.classes_ ) ] )

## 1.3.2 Eliminate duplicates by converting them into ones
pre2010_adj = pre2010_adj.tocsr( )
pre2010_adj.data = np.ones_like( pre2010_adj.data )

## 1.4.1 Find the vertices of the pre 2010 period that are in post-2010
post2010 = dblp[ dblp.year > 2010 ]
common_vertices = np.intersect1d( pre2010_values,
	np.union1d( post2010[ "author1" ].values, post2010[ "author2" ].values ) )

## 1.4.2 Remove completely new vertices from post 2010 data
post2010 = post2010[ (
    match( common_vertices, post2010[ "author1" ].values ) &
    match( common_vertices, post2010[ "author2" ].values ) ) ]
del common_vertices

## 1.4.3 Map the post 2010 vertices to pre 2010 vertices and construct the adjacency matrix.
for col in [ 'author1', 'author2', ] :
    post2010[ col ] = le.transform( post2010[ col ] )

## 1.4.4 The adjacency matrix
post2010_adj = sp.coo_matrix( ( np.ones( post2010.shape[0], dtype = np.bool ),
        ( post2010[ "author1" ].values, post2010[ "author2" ].values )
    ), shape = pre2010_adj.shape ).tolil( )

## 1.4.5 Leave only those edges in the post 2010 dataset, which had not existed during 2000-2010.
post2010_adj[ pre2010_adj.nonzero( ) ] = 0

## 1.4.6 Eliminate duplicate edges and transform into a CSR format
post2010_adj = post2010_adj.tocsr( )
post2010_adj.data = np.ones_like( post2010_adj.data )

print post2010_adj.__repr__( )
print pre2010_adj.__repr__( )

## 1.5.1 All edes of the post2010 graph are included and considered to be positive examples
positive = np.append( *( c.reshape((-1, 1)) for c in post2010_adj.nonzero( ) ), axis = 1 )

## 1.5.2 Generate a sample of vertex pairs with no edge in both periods
negative = np.random.choice( pre2010_adj.shape[ 0 ], size = ( 2 * positive.shape[0], positive.shape[1] ) )

## 1.5.3 Compie the final training dataset.
E = np.vstack( ( positive, negative ) )
y = np.vstack( ( np.ones( ( positive.shape[ 0 ], 1 ), dtype = np.float ),
                np.zeros( ( negative.shape[ 0 ], 1 ), dtype = np.float ) ) )

## So, finally, we got ourselves a trainig set of edges with 2:1 negative-to-postive ratio.

def phi_degree( edges, A ) :
	deg = A.sum( axis = 1 ).astype( np.float )
	return np.append( deg[ edges[ :, 0 ] ], deg[ edges[ :, 1 ] ], axis = 1 )

def __sparse_sandwich( edges, A, W = None ) :
    AA = A.dot( A.T ) if W is None else A.dot( W ).dot( A.T )
    result = AA[ edges[:,0], edges[:,1] ]
    del AA ; gc.collect( 0 ) ; gc.collect( 1 ) ; gc.collect( 2 )
    return result.reshape(-1, 1)

def phi_adamic_adar( edges, A ) :
    inv_log_deg = 1.0 / np.log( A.sum( axis = 1 ).getA1( ) )
    inv_log_deg[ np.isinf( inv_log_deg ) ] = 0
    result = __sparse_sandwich( edges, A, sp.diags( inv_log_deg, 0 ) )
    del inv_log_deg ; gc.collect( 0 ) ; gc.collect( 1 ) ; gc.collect( 2 )
    return result

def phi_common_neighbours( edges, A ) :
    return __sparse_sandwich( edges, A )

def __sparse_pagerank( A, beta = 0.85, one = None, niter = 1000, rel_eps = 1e-6, verbose = True ) :
## Initialize the iterations
	one = one if one is not None else np.ones( ( 1, A.shape[ 0 ] ), dtype = np.float )
	one = sp.csr_matrix( one / one.sum( axis = 1 ) )
## Get the out-degree
	out = np.asarray( A.sum( axis = 1 ).getA1( ), dtype = np.float )
## Obtain the mask of dangling vertices
	dangling = np.where( out == 0.0 )[ 0 ]
## Correct the out-degree for sink nodes
	out[ dangling ] = 1.0
## Just one iteration: all dangling nodes add to the importance of all vertices.
	pi = np.full( ( one.shape[0], A.shape[0] ), 1.0 / A.shape[ 0 ], dtype = np.float )
## If there are no dangling vertices then use simple iterations
	kiter, status = 0, -1
## Make a stochastic matrix
	P = sp.diags( 1.0 / out, 0, dtype = np.float ).dot( A ).tocsc( )
	while kiter < niter :
## make a copy of hte current ranking estimates
		pi_last = pi.copy( )
## Use sparse inplace operations for speed. Firstt the random walk part
		pi *= beta ; pi *= P
## Now the teleportaiton ...
		pi += ( 1 - beta ) * one
##  ... and dangling vertices part
		if len( dangling ) > 0 :
			pi += beta * one.multiply( np.sum( pi_last[ :, dangling ], axis = 1 ).reshape( ( -1, 1 ) ) )
## Normalize
		pi /= np.sum( pi, axis = 1 )
		if np.sum( np.abs( pi - pi_last ) ) <= one.shape[0] * rel_eps * np.sum( np.abs( pi_last ) ) :
			status = 0
			break
## Next iteration
		kiter += 1
		if kiter % 10 == 0 :
			print kiter
	return pi, status, kiter

def phi_gpr( edges, A, verbose = True ) :
	pi, s, k = __sparse_pagerank( A, one = None, verbose = verbose )
	return np.concatenate( ( pi[ :, edges[ :, 0 ] ], pi[ :,edges[ :, 1 ] ] ), axis = 0 ).T

def phi_ppr( edges, A ) :
	result = np.empty( edges.shape, dtype = np.float )

    return __sparse_sandwich( edges, A )

phi_12 = phi_degree( E, pre2010_adj )
phi_34 = phi_gpr( E, pre2010_adj )
phi_5 = phi_adamic_adar( E, pre2010_adj )
phi_6 = phi_common_neighbours( E, pre2010_adj )

X = np.hstack( ( phi_12, phi_34, phi_5, phi_6 ) )

from sklearn.cross_validation import train_test_split
X_modelling, X_main, y_modelling, y_main = train_test_split( X, y.ravel( ), train_size = 0.20 )

from sklearn.linear_model import LogisticRegression
clf = LogisticRegression(  ).fit( X_modelling, y_modelling )

theta = clf.predict_proba( X_main )

A = np.tensordot( np.vstack( ( y_main, 1-y_main ) ).T, np.log( np.clip( theta, 1e-15, 1-1e-15 ) ), ( 0, 0 ) )
logloss = -np.sum( np.diag( A ) / y_main.shape[ 0 ] )

from sklearn.cross_validation import cross_val_score
cross_val_score( clf, X_main, y_main, n_jobs = -1, verbose = 50, cv = 10 )

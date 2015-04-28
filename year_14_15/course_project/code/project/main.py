# -*- coding: UTF-8 -*-
from datetime import datetime

import numpy as np
from montecarlo import mc_run as montecarlo

from fgn import fbm
from crossing_tree import xtree_build

def mc_kernel( generate_sample, **op ) :
## Generate a replication of the process -- a sample path
	T, X = generate_sample( )
## G. Decrouez 2015-02-12: the selection of the spacing of the finest grid based on
##  the scale of the process is crucial as it allows comparison of the crossing tree
##  between different sample paths.
## 2015-04-08 : I think 6-level deep crossing tree is too poor for any analysis which
##  is why it is necessary to take the variance, i.e the inherent scale, into account.
	delta = np.std( np.diff( X ) )
## max_levels -- the number of level of the tree to construct (from the finest grid)
##  beyond this level (L+1, L+2, ...) the levels are pooled.
## max_crossing -- the threshold for the number of subcrossings beyond which they
##  are considered to be in the tail. The distribution data would contain bins:
##  {2}, {4}, ... {2k}, ... {2K, 2(K+1), ...}, where K = max_crossing
	max_levels, max_crossings = op.get( 'L', 6 ), op.get( 'K', 20 )
## Tnk[n][k] -- the time (approximate) of the k-th crossing of the gird with resolution
##  \delta 2^n. XTnk[n][k] -- the value of the process at the time of the crossing:
##  alwayts equal to the shifted and scaled grid level crossed.
## Znk[n][k] -- the number of subcrossings of a finer grid (crossings of
##  \delta 2^{n-1}) that make up the k-th corssing of a coarser grid (\delta 2^n).
##  Undefined (empty array) for n = 0.
## Vnk[n][k] -- the number of up-down and down-up excursions in the k-th
##  crossing of size \delta 2^n. Undefined for n = 0.
## Wnk[n][k] -- the waiting time between the k-th and k+1-st crossing of the grid
##  of spacing \delta 2^n
	Tnk, XTnk, Znk, Vnk, Wnk = xtree_build( T, X, delta = delta )
## Get the total number of crosssings of \delta 2^n resolution
## Nn[n] -- the total number of crossings of grid with spacing \delta 2^n
	Nn = np.array( [ len( Tk ) - 1 for Tk in Tnk ],
		dtype = np.float ).reshape( ( len( Tnk ), -1 ) )
## Dnm[n][m] -- m<M : the total number of crossings of grid \delta 2^{n+1} with 2(m+1)
##  subcrossings of grid \delta 2^n. The values in column M are the number
##  of crossings of \delta 2^{n+1} with not less than 2(M+1) subcrossings.
## 32bit Integers should be enough
	Dnm = np.zeros( ( max_levels + 1, max_crossings // 2 ), dtype = np.float )
	for n, Zk in enumerate( Znk[ 1: ], 0 ) :
		n = max_levels if n > max_levels else n
		Z_count, Z_freq = np.unique( Zk, return_counts = True )
## Truncate the observed crossings
		Z_count = np.minimum( Z_count, max_crossings )
		mask = ( Z_count < max_crossings )
		Dnm[ n, Z_count[ mask ] // 2 - 1 ] += Z_freq[ mask ]
		Dnm[ n, max_crossings // 2 - 1 ] += np.sum( Z_freq[ ~mask ] )
## VDn[n][d][e] -- the total number of up-down(e=0) and down-up(e=1) excursions
##  (/\ and \/ subcrosssings od \delta 2^n respectively) in a downward(d=0) or
##  upward(d=1) crossing of spacing \delta 2^{n+1} (level n+1). Levels beyond
##  max_levels are agregated.
	VDn = np.zeros( ( max_levels + 1, 2, 2 ), dtype = np.int )
	for n, Vk in enumerate( Vnk[ 1: ], 0 ) :
		n = max_levels if n > max_levels else n
		VDn[ n, 0 ] += np.sum( Vk[ Vk[ :, 2 ] < 0 ], axis = 0 )[:2]
		VDn[ n, 1 ] += np.sum( Vk[ Vk[ :, 2 ] > 0 ], axis = 0 )[:2]
	return Nn, Dnm, VDn


if __name__ == '__main__' :
	N = 2**18+1 ; M = 10
	for H in  np.linspace( .5, .95, num = 10 ) :
		P = int( np.log2( N - 1 ) )
		print "Monte carlo (%d) for FBM(2**%d+1, %.4f):" % ( M, P, H )
## Get the current timestamp
		run_dttm = datetime.utcnow( )
## Iinitalize the generator
		generator = fbm( N = N, H = H )
## Run the experiment
		result = montecarlo( generator, mc_kernel,
			replications = M, parallel = P < 21, debug = False,
			quiet = False, # processes = 2,
			L = 10, K = 20 )
## Create a meaningful name for the output data blob
		np.savez_compressed( "./output/fbm%s_%d_%.2f_%d" % (
				run_dttm.strftime( "%Y%m%d-%H%M%S" ), P, H, M ),
			Njn  = np.array( [ n for w, j, ( n, _, _ ) in result ] ),
			Djnm = np.array( [ d for w, j, ( _, d, _ ) in result ] ),
			VDjn = np.array( [ v for w, j, ( _, _, v ) in result ] ) )
## To access use: dat = np.load(..) ; dat['Djnm'], dat['Njn'], dat['VDjn']

## Generate default name for each Monte Carlo result
	# np.savez_compressed( './test.out', **{ "rep%04d"%(i,): d
	#  	for _, i, (d, _) in result } )

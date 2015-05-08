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
## 2015-04-29 : I think lower (finer) levels of the tree are too biased for the analysis
##  the bias is due to the discretization of the sample path.
	delta_method = op.get( 'delta', 'std' ).lower( )
	if delta_method == 'std' :
		delta = np.std( np.diff( X ) )
	elif delta_method == 'iqr' :
## Use the interquartile range
		delta = np.subtract( *np.percentile( np.diff( X ), [ 75, 25 ] ) )
	else :
		delta = 1.0
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
	Nn = np.zeros( ( max_levels + 1 + 1, 1 ), dtype = np.int )
	for n, Tk in enumerate( Tnk, 0 ) :
		n = max_levels + 1 if n > max_levels + 1 else n
## The correct number of crossings of level \delta 2^n is the number of consecutive
##  pairs with Tnk_i < Tnk_{i+1}.
		Nn[ n ] += len( Tk ) - 1
## Dnk[n][k] -- k<K : the total number of crossings of grid \delta 2^{n+1} with 2(k+1)
##  subcrossings of grid \delta 2^n. The values in column K are the number
##  of crossings of \delta 2^{n+1} with not less than 2(K+1) subcrossings.
## 32bit Integers should be enough
	Dnk = np.zeros( ( max_levels + 1, max_crossings // 2 ), dtype = np.int )
	for n, Zk in enumerate( Znk[ 1: ], 0 ) :
		n = max_levels if n > max_levels else n
		Z_count, Z_freq = np.unique( Zk, return_counts = True )
## Truncate the observed crossings
		Z_count = np.minimum( Z_count, max_crossings )
		mask = ( Z_count < max_crossings )
		Dnk[ n, Z_count[ mask ] // 2 - 1 ] += Z_freq[ mask ]
		Dnk[ n, max_crossings // 2 - 1 ] += np.sum( Z_freq[ ~mask ] )
## Vnde[n][d][e] -- the total number of up-down(e=0) and down-up(e=1) excursions
##  (/\ and \/ subcrosssings od \delta 2^n respectively) in a downward(d=0) or
##  upward(d=1) crossing of spacing \delta 2^{n+1} (level n+1). Levels beyond
##  max_levels are agregated.
	Vnde = np.zeros( ( max_levels + 1, 2, 2 ), dtype = np.int )
	for n, Vk in enumerate( Vnk[ 1: ], 0 ) :
		n = max_levels if n > max_levels else n
		Vnde[ n, 0 ] += np.sum( Vk[ Vk[ :, 2 ] < 0 ], axis = 0 )[:2]
		Vnde[ n, 1 ] += np.sum( Vk[ Vk[ :, 2 ] > 0 ], axis = 0 )[:2]
	return Nn, Dnk, Vnde

if __name__ == '__main__' :
	N = 2**23+1 ; M = 2000
	for delta_method in [ 'iqr', 'std' ] :
		for H in np.linspace( .5, .95, num = 10 ) :
			P = int( np.log2( N - 1 ) )
			print "Monte carlo (%d) for FBM(2**%d+1, %.4f):" % ( M, P, H )
## Get the current timestamp
			run_dttm = datetime.utcnow( )
## Iinitalize the generator
			generator = fbm( N = N, H = H )
## Run the experiment
			result = montecarlo( generator, mc_kernel,
				processes = 4, debug = False, quiet = False, parallel = True,
				replications = M, delta = delta_method, L = 20, K = 40 )
## Create a meaningful name for the output data blob
			np.savez_compressed( "C:/Users/ivannz/Dropbox/study_notes/year_14_15/course_project/code/output/fbm_%s_%s_%d_%.4f_%d" % (
					delta_method.lower( ), run_dttm.strftime( "%Y%m%d-%H%M%S" ), P, H, M ),
				Njn   = np.array( [ n for w, j, ( n, _, _ ) in result ] ),
				Djnk  = np.array( [ d for w, j, ( _, d, _ ) in result ] ),
				Vjnde = np.array( [ v for w, j, ( _, _, v ) in result ] ) )
## To access use: dat = np.load(..) ; dat['Djnk'], dat['Njn'], dat['Vjnde']
## For analysis:
##  Vnd = np.sum( Vnde, axis = 2, dtype = np.float ).reshape( VDn.shape[:2] + ( 1, ) )
##  Vnde / Vnd <- the conditional distribution of excursions
##  Dn = np.sum( Dnm, axis = 1, dtype = np.float ).reshape( Dnm.shape[:1] + ( 1, ) )
##  Dnk / Dn <- the subcrossing number distribution
## Generate default name for each Monte Carlo result
	# np.savez_compressed( './test.out', **{ "rep%04d"%(i,): d for _, i, (d, _) in result } )

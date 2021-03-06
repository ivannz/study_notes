# -*- coding: UTF-8 -*-
## Base modules
import os, re
import numpy as np
from datetime import datetime

## Project modules
from fgn import fbm
from hermite import hermite
from weierstrass import weierstrass
from montecarlo import mc_run as montecarlo
from crossing_tree import xtree_build

## This procedure runs a single replication of the Monte-Carlo experiment.
def mc_kernel( replication, generate_sample, **op ) :
## Generate a replication of the process -- a sample path
	T, X = generate_sample( )
	return path_kernel( T, X, **op )

## This routine just unpacks the arguments and performs necessary pre computations
def path_kernel( T, X, **op ) :
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
	elif delta_method == 'rng' :
## Use the range estimate as suggested by Geoffrey on 2015-05-28
		delta = np.subtract( *np.percentile( X, [ 100, 0 ] ) ) * 2**( -6 )
	elif delta_method == 'med' :
## Use the median as suggested in [Jones, Rolls; 2009] p. 11 (0911.5204v2)
		delta = np.median( np.abs( np.diff( X ) ) )
	else :
## By default, delta is set to 1, which is the most inferior choice, since it disregards
##  the base scale of the sample path.
		delta = 1.0
	max_levels, max_crossings = op.get( 'L', 6 ), op.get( 'K', 20 )
## Analyse the sample path
	return path_analyse( T, X, delta, max_levels, max_crossings )

## This procedure runs a single replication of the Monte-Carlo experiment.
## delta -- the base scale of the crossing tree. The larger, the more biased the lower
##  tree levels are.
## max_levels -- the number of levels of the tree to construct (from the finest grid up)
##  Beyond this level (L+1, L+2, ...) the levels are pooled.
## max_crossings -- the threshold for the number of subcrossings beyond which they are
##  considered to be in the "tail". The distribution data would contain bins :
##  [ {2}, {4}, ... {2k}, ... {2K, 2(K+1), ...}], where K = max_crossing
def path_analyse( T, X, delta = 1.0, max_levels = 6, max_crossings = 20 ) :
## Tnk[n][k] -- the time (approximate) of the k-th crossing of the gird with resolution
##  \delta 2^n. Xnk[n][k] -- the value of the process at the time of the crossing:
##  alwayts equal to the shifted and scaled grid level crossed.
## Znk[n][k] -- the number of subcrossings of a finer grid (crossings of
##  \delta 2^{n-1}) that make up the k-th crossing of a coarser grid (\delta 2^n).
##  Undefined (empty array) for n = 0.
## Vnk[n][k] -- the number of up-down and down-up excursions in the k-th crossing
##  of size \delta 2^n. Undefined for n = 0. The format is (# of /\, # of \/, ±1).
## Wnk[n][k] -- the waiting time between the k-th and k+1-st crossing of the grid
##  of spacing \delta 2^n
	Tnk, Xnk, Znk, Vnk, Wnk = xtree_build( T, X, delta = delta )
## Get the total number of crossings of \delta 2^n resolution
## Nn[n] -- the total number of crossings of grid with spacing \delta 2^n
	Nn = np.zeros( ( 1 + max_levels + 1, 1 ), dtype = np.int )
	for n, Xk in enumerate( Xnk, 0 ) :
		n = max_levels + 1 if n > max_levels + 1 else n
## The correct number of crossings of level \delta 2^n is the number of consecutive
##  pairs with Tnk_i < Tnk_{i+1}.
		Nn[ n ] += len( Xk ) - 1
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
##  max_levels are aggregated.
	Vnde = np.zeros( ( max_levels + 1, 2, 2 ), dtype = np.int )
	for n, Vk in enumerate( Vnk[ 1: ], 0 ) :
		n = max_levels if n > max_levels else n
		Vnde[ n, 0 ] += np.sum( Vk[ Vk[ :, 2 ] < 0 ], axis = 0 )[:2]
		Vnde[ n, 1 ] += np.sum( Vk[ Vk[ :, 2 ] > 0 ], axis = 0 )[:2]
## Make a crude summary of crossing durations: 
	prc = np.array( [ 0.5, 1.0, 2.5, 5.0, 10, 25, 50, 75, 90, 95, 97.5, 99, 99.5 ] )
## Wnp[n][p] -- the p-th empirical quantile of the n-th level crossing durations.
##  Crossing durations are not aggregated in the last row of the output.
	Wnp = np.zeros( ( max_levels, ) + prc.shape, dtype = np.float )
## The average crossing duration and its standard deviation
	Wbarn = np.zeros( ( max_levels, 1 ), dtype = np.float )
## The average crossing duration and its standard deviation
	Wstdn = np.zeros( ( max_levels, 1 ), dtype = np.float )
	for n, Wk in enumerate( Wnk[1:], 0 ) :
		if len( Wk ) and n < max_levels :
## Get the mean, standard deviation and the quantiles of non-empty levels only.
			Wbarn[ n ], Wstdn[ n ], Wnp[ n ] = np.average( Wk ), np.std( Wk ), np.percentile( Wk, prc )
	return Nn, Dnk, Vnde, ( Wnp, Wbarn, Wstdn )

## A procedure for loading results of a saved simulation
def sim_load( file_name, return_durations = False ) :
## Load the npz file
	data = np.load( file_name )
## Remove the extension and split by underscores: the format is (G, D, dttm, P, H, M)
	par = re.sub( r'(.*)\.npz$', r'\1', os.path.basename( file_name ) ).split( '_' )
	if not return_durations :
## return the basic results of the simulation : tuple( H, Njn, Djnk, Vjnde )
		return float( par[ -2 ] ), data[ 'Njn' ], data[ 'Djnk' ], data[ 'Vjnde' ]
	else :
## Get the dimensions of the data
		M, L, K = data[ 'Djnk' ].shape
## Prepare the data on crossing durations
		Wjnp   = data[ 'Wjnp' ]   if 'Wjnp'   in data else np.empty( ( M, L, 0 ) )
		Wbarjn = data[ 'Wbarjn' ] if 'Wbarjn' in data else np.empty( ( M, L, 0 ) )
		Wstdjn = data[ 'Wstdjn' ] if 'Wstdjn' in data else np.empty( ( M, L, 0 ) )
## return the results of the simulation : basic plus summary of duration data
		return float( par[ -2 ] ), data[ 'Njn' ], data[ 'Djnk' ], data[ 'Vjnde' ], Wjnp, Wbarjn, Wstdjn

## A procedure for saving the results of the Monte-Carlo simulation
def sim_save( filename, result, save_durations = False, **extra ) :
	npz_data = dict( **extra )
	npz_data[ "Njn" ]     = np.array( [ n for wrk, j, ( n, _, _, _ ) in result ] )
	npz_data[ "Djnk" ]    = np.array( [ d for wrk, j, ( _, d, _, _ ) in result ] )
	npz_data[ "Vjnde" ]   = np.array( [ v for wrk, j, ( _, _, v, _ ) in result ] )
## If requested save the measured parameters of durations as well
	if save_durations :
		npz_data["Wjnp"]    = np.array( [ w for wrk, j, ( _, _, _, ( w, _, _ ) ) in result ] )
		npz_data["Wbarjn"]  = np.array( [ b for wrk, j, ( _, _, _, ( _, b, _ ) ) in result ] )
		npz_data["Wstdjn"]  = np.array( [ s for wrk, j, ( _, _, _, ( _, _, s ) ) in result ] )
	np.savez_compressed( filename, **npz_data )

## A service procedure for dumping a directory
def list_files( path = './', pattern = r'\.npz$' ) :
## If the regular expression fails to compile, return an empty list of files,
##  since no file name would match an invalid expression.
	try :
		mask = re.compile( pattern, re.UNICODE )
## Return full path to located files
		return [ os.path.join( os.path.realpath( path ), filename )
			for filename in os.listdir( path ) if mask.search( filename ) ]
	except :
		return [ ]

if __name__ == '__main__' :
## The base path of the the results
	basepath = os.path.realpath( "./output/final" )
####################################################################################################
	if False :
## Fractional Brownian Motion
		N, M = 2**18+1, 1000
		P = int( np.log2( N - 1 ) )
## Loop over the base scale methods
		for delta_method in [ 'med', ] :
## Create the necessary directory structure
			target_path = os.path.join( basepath, "FBM_%d" % ( P, ), delta_method )
			os.makedirs( target_path )
## Run simulations for different Hurst exponents
			for H in np.linspace( .5, .9, num = 5 ) :
## Initalize the generator
				generator = fbm( N = N, H = H, time = True )
				print "Monte-Carlo (%d) for FBM(2**%d+1, %.4f), %s:" % ( M, P, H, delta_method, )
## Get the current timestamp
				start_dttm = datetime.utcnow( )
## Run the experiment
				results = montecarlo( generator, mc_kernel,
					processes = 2, debug = False, quiet = False, parallel = True,
					replications = M, delta = delta_method, L = 20, K = 40 )
## Create a meaningful name for the output data blob
				file_name = "FBM_%s_%s_%d_%.4f_%d" % ( delta_method.lower( ),
					start_dttm.strftime( "%Y%m%d-%H%M%S" ), P, H, M )
## Save the data blob
				sim_save( os.path.join( target_path, file_name ), results, save_durations = True )
####################################################################################################
	if False :
## Hermite processes
## The parameters of the simulation
		N, K, M = 2**18+1, 2**4, 1000
		P = int( np.log2( N - 1 ) )
		for delta_method in [ 'med', ] :
			for D in [ 2, 3, 4, ] :
				target_path = os.path.join( basepath,
					"HRM-%d_%d-%d" % ( D, P, K, ), delta_method )
				os.makedirs( target_path )
				for H in np.linspace( .6, .9, num = 4 ) :
## Initialize the Hermite process generator
					generator = hermite( N = N, d = D, H = H, K = K, time = True )
					print "Monte-Carlo (%d) for HRM-%d(2**%d+1-%d, %.4f), %s:" % ( M, D, P, K, H, delta_method, )
					start_dttm = datetime.utcnow( )
					results = montecarlo( generator, mc_kernel,
						processes = 2, debug = False, quiet = False, parallel = True,
						replications = M, delta = delta_method, L = 20, K = 40 )
					file_name = "HRM-%d_%s_%s_%d-%d_%.4f_%d" % ( D, delta_method.lower( ),
						start_dttm.strftime( "%Y%m%d-%H%M%S" ), P, K, H, M )
					sim_save( os.path.join( target_path, file_name ), results, save_durations = True )
####################################################################################################
	if True :
## Weierstrass processes
		N, M = 2**18+1, 1000
		P = int( np.log2( N - 1 ) )
## Loop over the base scale methods
		for delta_method in [ 'med', 'iqr', 'rng', ] :
## Create the necessary directory structure
			target_path = os.path.join( basepath, "WEI_%d" % ( P, ), delta_method )
			os.makedirs( target_path )
## Run simulations for different Hurst exponents
			for H in [0.5] :#np.linspace( .5, .9, num = 5 ) :
## Initalize the generator
				generator = weierstrass( N = N, H = H, nu0 = 1.2, nue = 1000 )
				print "Monte-Carlo (%d) for WEI(2**%d+1, %.4f), %s:" % ( M, P, H, delta_method, )
## Get the current timestamp
				start_dttm = datetime.utcnow( )
## Run the experiment
				results = montecarlo( generator, mc_kernel,
					processes = 2, debug = False, quiet = False, parallel = True,
					replications = M, delta = delta_method, L = 20, K = 40 )
## Create a meaningful name for the output data blob
				file_name = "WEI_%s_%s_%d_%.4f_%d" % ( delta_method.lower( ),
					start_dttm.strftime( "%Y%m%d-%H%M%S" ), P, H, M )
## Save the data blob
				sim_save( os.path.join( target_path, file_name ), results, save_durations = True )






####################################################################################################
####################################################################################################
## To access use: dat = np.load(..) ; dat['Djnk'], dat['Njn'], dat['Vjnde']
## For analysis:
##  Vnd = np.sum( Vnde, axis = 2, dtype = np.float ).reshape( VDn.shape[:2] + ( 1, ) )
##  Vnde / Vnd <- the conditional distribution of excursions
##  Dn = np.sum( Dnm, axis = 1, dtype = np.float ).reshape( Dnm.shape[:1] + ( 1, ) )
##  Dnk / Dn <- the subcrossing number distribution
## Generate default name for each Monte Carlo result
## 	np.savez_compressed( './test.out', **{ "rep%04d"%(i,): d for _, i, (d, _) in result } )

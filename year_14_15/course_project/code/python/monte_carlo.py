#! /usr/bin/env python
# -*- coding: UTF-8 -*-

import numpy as np

from crossing_tree import xtree_build

## The Monte Carlo kernel performs a single experiment and returns its
## K'xT slice of the result.
def monte_carlo_kernel( generator, K, T ) :
## K is the maximal number of subcrossings beyond which the data
##  on subcrossings is agregateed into the tail. (truncation parameter)
	K = 2 * ( K // 2 + 1 )
## T is the number of detailed levels of the crossing tree to stored in
##  output. Any crossings of grids coarser than the treshold are
##  aggregated.
	distr = np.zeros( ( K // 2, T + 1 ), dtype = np.int )
## Generate a replication of the process
	t, x = generator( )
## Get the hitting times and poitns
	ht, hp, subx = xtree_build( t, x )
## Run the chi_square test on each level of the tree
##	ct = ctable( subx[ :-1 ], subx[ 1: ] )
##	chisq_test( ct )
	for level, xing in enumerate( subx[1:], 0 ) :
## Count the absolute frequencies
		c, f = np.unique( xing, return_counts = True )
## We collect the data on the distribution truncated by K (including).
## Thus truncate the number of subcrossings by K -- everything less is
##  detailed, otherwise -- agregated into tail
		c = np.minimum( c, K )
## Similarly truncate the level by T + 1
		if level >= T : level = T
## Fill the offspring distribution matrix
		for c, f in zip( c, f ) :
			distr[ c // 2 - 1, level ] += f
## Use straightforward bivariate regression implementation.
## We could also use ML estimators, which is much better!
## See the handwritten notes.
	return distr

def monte_carlo_serial( generator, M = 100, K = 4, T = 2 ) :
## K is the maximal number of subcrossings beyond which the data
##  on subcrossings is agregateed into the tail. (truncation parameter)
	K = 2 * ( K // 2 + 1 )
## T is the number of detailed levels of the crossing tree to stored in
##  output. Any crossings of grids coarser than the treshold are
##  aggregated.
	distr = np.zeros( ( K // 2, T + 1, M ), dtype = np.int32 )
## The output of the ChSquared independence test
# 	chisq = np.empty( ( ) )
## The parallel computing part: use synchronous computations
#	cli = mp.Client() ; cli.block = True
#	clu = clients.load_balanced_view()
#	clu.map(monte_karlo_kernel, [5, 6, 7, 8], [8, 9, 10, 11]        
	for m in xrange( M ) :
		if m % ( M // 10 ) == 0 : print "%d " % ( m, )
## Generate a replication of the process
		t, x = generator( )
## Get the hitting times and poitns
		ht, hp, subx = xtree_build( t, x )
## Run the chi_square test on each level of the tree
# 		ct = ctable( subx[ :-1 ], subx[ 1: ] )
##		chisq_test( ct )
## Pool together the second and the third layer subcrossings
		pool = []
		for level, xing in enumerate( subx[1:], 0 ) :
## Count the absolute frequencies
			c, f = np.unique( xing, return_counts = True )
## We collect the data on the distribution truncated by K (including).
## Thus truncate the number of subcrossings by K -- everything less is
##  detailed, otherwise -- agregated into tail
			c = np.minimum( c, K )
## Similarly truncate the level by T + 1
			if level >= T : level = T
## Fill the offspring distribution matrix
			for c, f in zip( c, f ) :
				distr[ c // 2 - 1, level, m ] += f
## Distr is KxTxM
	return distr

## DD = monte_carlo_serial( genr, M=100, K = 16, T = 3 )
def monte_carlo_parallel( genr, M = 100, K = 16, T = 3 ) :
    import time
## Setup the cluster for the monte carlo experiment.
    cluster = mp_mc_setup( genr, M, K, T )
## Run the expriment asynchronously
    result = cluster.apply_async( mp_mc_worker )
## Wait intil the jobs are complete
    while not result.ready( ):
        time.sleep( 2 )
## Track the porgress
        print result.progress
    return np.concatenate([r for r in result if r])
## Initialize the cluster to run the Monte Carlo experiment
import IPython.parallel as mp
def mp_mc_setup( genr, M = 100, K = 4, T = 2 ) :
## The parallel computing part: use synchronous computations
	cli = mp.Client()
	cli.clear()
	clu = cli.direct_view()
## Make the workers import the necessary dependencies
	clu.execute('import numpy as np', block = True )
## The crossing tree toolkit
	clu.execute('from crossing_tree import xtree_build', block = True )
# 	with clu.sync_imports( ) :
# 		import numpy
# 		from monte_carlo import monte_carlo_kernel
# 		from crossing_tree import xtree_build
## The generators
# 		from Weierstrass import synth_Weier
# 		from synthfbmcircul import synth_fbm, synth_fgn
# 		from hsssi_processes import synth_Rosenblatt, synth_Hermite3, synth_Hermite4
# Distribute the workload evenly among the members of the cluster
	clu.scatter( 'local_replications', range( M ) )
## Pass the necesary environment to the cluster.
	clu.push({
		'generator' : genr, 'K' : K, 'T': T,
		'M': M, 'kernel': monte_carlo_kernel })
	return clu

## This is a genral procedure to be run on each node of the cluster.
def mp_mc_worker( ) :
# 	return ( min(local_replications), max(local_replications) )
## Define a local storage for the results of calls to the kernel
#     result = []
#     for m in local_replications :
#         distr = kernel( m, genr, K, T )
#         result.append( distr )
## Perform the appointed number of monte carlo local_replications
	return [ kernel( generator, K, T ) for m in local_replications ]
#! /usr/bin/env python
# -*- coding: UTF-8 -*-

import numpy as np

# from crossing_tree import xtree_build

import IPython.parallel as mp
from IPython.parallel import interactive

import struct as struct

def monte_carlo_serial( generator, kernel, M = 100, quiet = False, **kwargs ) :
## Return a MxKxT 3D array made from M concatenated 2D KxT slices
	import time
	tic = time.time( )
	result = [ kernel( generator, **kwargs ) for m in xrange( M ) ]
	if quiet :
		print( "--: %.3f" % ( time.time( ) - tic ) )
	return np.concatenate( [ result ] )

## DD = monte_carlo_serial( generator, M=100, K = 16, T = 3 )
def monte_carlo_parallel( generator, kernel, M = 100, quiet = False, **kwargs ) :
	import time
## Setup the cluster for the monte carlo experiment.
	cluster = __mp_mc_setup( generator, kernel, M, **kwargs )
	if quiet :
## Run the expriment synchronously
		result = cluster.apply_sync( __mp_mc_worker )
	else :
## Run the expriment asynchronously
		result = cluster.apply_async( __mp_mc_worker )
## Wait intil the jobs are complete
		tic = time.time( )
		while not result.ready( ) :
			time.sleep( 1 )
## Track the porgress
			print( "%i: %.3f" % ( result.progress, time.time( ) - tic ) )
	cluster.clear( block = True )
	return np.concatenate([r for r in result if r])

## This is a genral procedure to be run on each node of the cluster.
@interactive
def __mp_mc_worker( ) :
## It is necessary to reininialize the random number generator so
##  as to avoid rndom number sequence collisions on different worker
##  units. cf. http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=5547156
##  http://stackoverflow.com/questions/12915177/same-output-in-different-workers-in-multiprocessing
	np.random.seed( seed = seed )
## Each process generator must be equipped with a reset method() which
##  resets its internal state.
	generator.reset( )
## Define a local storage for the results of calls to the kernel
	return [ kernel( generator, **local_kwargs ) for m in local_replications ]
## http://stackoverflow.com/questions/10857250/python-name-space-issues-with-ipython-parallel/10859394#10859394

## Initialize the cluster to run the Monte Carlo experiment
def __mp_mc_setup( generator, kernel, M = 100, **kwargs ) :
## The parallel computing part: use synchronous computations
	cli = mp.Client( )
	clu = cli.direct_view( )
	clu.clear( block = True )
## Make the workers import the necessary dependencies
	clu.execute( 'import numpy as np', block = True )
## The crossing tree toolkit
	clu.execute( 'from crossing_tree import xtree_build', block = True )
## The generators
	clu.execute( 'from Weierstrass import synth_Weier', block = True )
	clu.execute( 'from synthfbmcircul import synth_fbm, synth_fgn', block = True )
	clu.execute( 'from hsssi_processes import synth_Rosenblatt, synth_Hermite3, synth_Hermite4', block = True )
# Distribute the workload evenly among the members of the cluster
	clu.scatter( 'local_replications', range( M ), block = True )
## The main problem is that it is possible that each worker reads
##  its seed from the same source and at the same time, which would
##  produces dangerously correlated results, even meaningless. Therefore
##  before starting the kernels, let the parent process generate
##  some entropy for each child process.
##  http://stackoverflow.com/questions/2396209/best-seed-for-parallel-process
## Generate 32bit values uniformly at random.
	seeds = np.random.randint( 0xFFFFFFFF, size = len( cli ) )
## Alternative way for UNIX systems is to read from /dev/random.
##	with open( "/dev/random", "rb" ) as dev :
##		seeds = [ struct.unpack( 'I', dev.read( 4 ) )[ 0 ] for c in cli ]
## Dispatch individual seed values to each worker.
	for c, seed in zip( cli, seeds ) :
		c.push( { 'seed' : seed } )
## Pass the necesary environment to the cluster.
	clu.push({ 'generator' : generator, 'local_kwargs' : kwargs, 'local_M': M,
		'kernel': interactive( kernel ) }, block = True )
	return clu

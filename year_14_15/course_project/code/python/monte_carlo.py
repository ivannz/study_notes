#! /usr/bin/env python
# -*- coding: UTF-8 -*-

import numpy as np

# from crossing_tree import xtree_build

import IPython.parallel as mp
from IPython.parallel import interactive

def monte_carlo_serial( generator, kernel, M = 100, **kwargs ) :
## Return a MxKxT 3D array made from M concatenated 2D KxT slices
	import time
	tic = time.time( )
	result = [ kernel( generator, **kwargs ) for m in xrange( M ) ]
	print( "--: %.3f" % ( time.time( ) - tic ) )
	return np.concatenate( [ result ] )

## DD = monte_carlo_serial( generator, M=100, K = 16, T = 3 )
def monte_carlo_parallel( generator, kernel, M = 100, **kwargs ) :
	import time
## Setup the cluster for the monte carlo experiment.
	cluster = __mp_mc_setup( generator, kernel, M, **kwargs )
## Run the expriment asynchronously
	result = cluster.apply_async( __mp_mc_worker )
## Wait intil the jobs are complete
	tic = time.time( )
	while not result.ready( ) :
		time.sleep( 1 )
## Track the porgress
		print( "%i: %.3f" % ( result.progress, time.time( ) - tic ) )
	return np.concatenate([r for r in result if r])

## This is a genral procedure to be run on each node of the cluster.
@interactive
def __mp_mc_worker( ) :
## Define a local storage for the results of calls to the kernel
	return [ kernel( generator, **local_kwargs ) for m in local_replications ]
## http://stackoverflow.com/questions/10857250/python-name-space-issues-with-ipython-parallel/10859394#10859394

## Initialize the cluster to run the Monte Carlo experiment
def __mp_mc_setup( generator, kernel, M = 100, **kwargs ) :
## The parallel computing part: use synchronous computations
	cli = mp.Client()
	cli.clear( )
	clu = cli.direct_view()
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
## Pass the necesary environment to the cluster.
	clu.push({ 'generator' : generator, 'local_kwargs' : kwargs, 'local_M': M,
		'kernel': interactive( kernel ) }, block = True )
	return clu

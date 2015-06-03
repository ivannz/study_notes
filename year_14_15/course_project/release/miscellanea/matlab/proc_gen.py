# -*- coding: UTF-8 -*-
## This should be run from the location of main.py
## Base modules
import os, gc
from datetime import datetime
import scipy.io as io
import numpy as np
from numpy.random import RandomState

## Project modules
from processes.fgn_pyfftw import fgn
from processes.fbm import fbm
from processes.hermite import hermite
from processes.weierstrass import weierstrass

def kernel( generator, **op ) :
	files = list( )
## Get the parameters
	path = os.path.realpath( op.get( 'path', './output' ) )
	prefix = op.get( 'prefix', 'default' )
	replications = op.get( 'replications', 10 )
	for m in xrange( replications ) :
## Generate a sample path
		T, X = generator( )
## Get the current timestamp
		gen_dttm = datetime.utcnow( )
## Make a sensible name for the file
		filename = "%s_%s_%04X.mat" % ( prefix,
			gen_dttm.strftime( "%Y%m%d-%H%M%S" ), np.random.randint( 65536 ), )
## Save into a matlab file
		io.savemat( os.path.join( path, filename ), { 'T': T, 'X': X } )
		del T, X
		gc.collect( )
	return files

if __name__ == '__main__' :
	path = os.path.realpath( './output' )
## Parameters of the process to be generated
	N, K, M = 2**16 + 1, 2**5, 20
## Iterate over the degrees and the hurst indices
	# for D in [ 2, 3, 4 ] :
	for D in [ 3, 4 ] :
		for H in np.linspace( .5, .9, num = 5 ) :
## Make the necessary directory structure
			target_path = os.path.join( path, "HRM_%d_%d-%d" % (
				D, int( np.log2( N - 1 ) ), K, ),
				"%.4f" % (H, ) )
			os.makedirs( target_path )
## Initialize the generator: process and randomness
			gen = hermite( fgn, N, D, H = H, K = K )
			gen.initialize( RandomState( ) )
## Generate the data
			kernel( gen, path = target_path, replications = M,
				prefix = "HRM_%d_%d-%d_%.4f" % ( D,
					int( np.log2( N - 1 ) ), K, H, ) )

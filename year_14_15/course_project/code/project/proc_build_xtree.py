# -*- coding: UTF-8 -*-
## This is a very poor and ill informed solution: disk IO had better
##  be sequential!

## Base modules
import os, re, gc
import numpy as np
import scipy.io as io
from datetime import datetime

## Project modules
from gen import empty
from main import path_kernel, sim_save, list_files
from montecarlo import mc_run as montecarlo

def offline_kernel( file, generator, **op ) :
	mat = io.loadmat( file )
	X, T = mat['X'][0], mat['T']
	del mat
	gc.collect( )
	if len( T ) < 1 :
		T = np.arange( len( X ), dtype = np.float ) / ( len( X ) - 1 )
	else :
		T = T[ 0 ]
	return path_kernel( T, X, **op )

if __name__ == '__main__' :
	# basepath = os.path.join( '.', 'output', 'HRM_2_20-32' )
	basepath = os.path.join( '.', 'output', 'HRM_3_20-32' )
	path = [ os.path.join( basepath, H ) for H in [
		'0.5000', '0.6000', '0.7000', '0.8000',  ] ]
	for p in path :
## Go!
		files = list_files( p, pattern = r"\.mat$" )
		print( "Sweeping through %s: %d files found." %( p, len( files ), ) )
## Pick the first file and extract as much information as possible from its filename
		prefix, degree, size, hurst = os.path.basename( files[0] ).split( "_" )[:4]
		prefix += "-%s"%( degree, )
## Loop thorugh all the possible base scale calculation methods
		for delta_method in [ 'std', 'iqr', 'med', ] :
## Get the current timestamp
			run_dttm = datetime.utcnow( )
## Run the analyzer in parallel
			result = montecarlo( empty( ), offline_kernel,
				processes = 1, debug = False, quiet = False, parallel = False,
				replications = files, delta = delta_method, L = 20, K = 40 )
## Get the datetime after the simulation has finished
			end_dttm = datetime.utcnow( )
## Create a meaningful name for the output data blob
			sim_save( os.path.join( basepath, "%s_%s_%s_%s_%s_%d" % ( prefix, delta_method.lower( ),
				run_dttm.strftime( "%Y%m%d-%H%M%S" ), size, hurst, len( files ) ) ),
				result, save_durations = True )

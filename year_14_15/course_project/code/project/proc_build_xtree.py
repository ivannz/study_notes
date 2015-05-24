# -*- coding: UTF-8 -*-
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
	delta_method, path = 'med', './output/raw/'
## Go!
	files = list_files( path, pattern = r"\.mat$" )
## Pick the first file and extract as much information as possible from its filename
	prefix, degree, size, hurst = os.path.basename( files[0] ).split( "_" )[:4]
	prefix += "-%s"%( degree, )
## Get the current timestamp
	run_dttm = datetime.utcnow( )
## Run the analyzer in parallel
	result = montecarlo( empty( ), offline_kernel,
		processes = 2, debug = False, quiet = False, parallel = True,
		replications = files, delta = delta_method, L = 20, K = 40 )
## Get the datetime after the simulation has finished
	end_dttm = datetime.utcnow( )
## Create a meaningful name for the output data blob
	sim_save( "./output/%s_%s_%s_%s_%s_%d" % ( prefix, delta_method.lower( ),
		run_dttm.strftime( "%Y%m%d-%H%M%S" ), size, hurst, len( files ) ),
		result, save_durations = True )

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
basepath = os.path.join( 'C:\\', 'TEMP', 'tp89957cff_f6cf_463c_b154_eac2aabe2c78', 'HRM-2_20-16' )

path = [ os.path.join( basepath, H ) for H in [
	'0.6', '0.7', '0.8', '0.9', ] ]

delta_method = 'med'

for p in path :
	files = list_files( p, pattern = r"\.mat$" )
## Pick the first file and pares its name
	print( "Sweeping through %s: %d files found." %( p, len( files ), ) )
## Pick the first file and extract as much information as possible from its filename
	prefix, size, hurst = os.path.basename( files[ 0 ] ).split( "_" )[:3]
	hurst = "%0.5f" % ( float( hurst ), )
## Get the current timestamp
	run_dttm = datetime.utcnow( )
	result = list( )
	for f in files :
## Load the Matlab data
		mat = io.loadmat( f )
		X = mat['data'][0]
		del mat
## Reconstruct the time
		T = np.arange( len( X ), dtype = np.float ) / ( len( X ) - 1 )
## Construct the crossing tree
		result.append( ( -1, -1, path_kernel( T, X, delta = delta_method, L = 20, K = 40 ) ) )
		gc.collect( )
## Save the analysis
	sim_save( os.path.join( basepath, "%s_%s_%s_%s_%s_%d" % ( prefix, delta_method.lower( ),
		run_dttm.strftime( "%Y%m%d-%H%M%S" ), size, hurst, len( files ) ) ),
		result, save_durations = True )
	del result[:]
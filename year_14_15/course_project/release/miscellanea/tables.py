# -*- coding: UTF-8 -*-
## Base modules
import re, os
import numpy as np

import warnings
warnings.filterwarnings( "ignore" )

import matplotlib.pyplot as plt

## Project modules
from main import sim_load, list_files
import self_similarity as ss

import plots as pp
base_path = os.path.realpath( "./output/release/10000-17" )
process_folders = [ 'FBM_17', 'HRM-2_17-16', 'HRM-3_17-16', 'HRM-4_17-16', 'WEI_17', ]

# base_path = os.path.realpath( "./output/release/1000-21" )
# process_folders = [ 'FBM_21', ]

## Quickly extract the info
kinds = [ process.split( '_' )[ 0 ] for process in process_folders ]

## Method name and the tree levels to use (chosen heuristically).
methods = [ ( 'med', 6, 7 ), ]
# methods = [ ( 'med', 6, 7 ), ]

results = dict( )
for method, _, _ in methods :
	for process in process_folders :
		path = os.path.join( base_path, process, method )
		if not os.path.exists( path ) :
			continue
## Extract the name of the process from its folder name
		kind = process.split( '_' )[ 0 ]
## Load the results results
		results_method = results.get( method, dict( { } ) )
		results_method[ kind ] = pp.load_files( list_files( path ) )
		results[ method ] = results_method


## Generate the body of a table for H = 0.6
def conjecture_table( method, hurst, body = "" ): 
	for kind in kinds :
		row = ""
## get the results
		simd = pp.find_hurst( method, kind, hurst, results )
		if simd is None :
			continue
		L, H, Njn, Djnk, Vjnde = simd
## Conjectured
		theta = 2**( 1 - 1 / H )
		row += "\multirow{2}{*}{%s-$%0.2f$}" %( kind, H, ) + " & $%0.3f$ & $%0.3f$ & $%0.3f$ & $%0.3f$ \\\\ \\cline{2-5} \n" % (
		theta, theta * (1-theta), theta * (1-theta)**2, theta * (1-theta)**3 )
## Empirical
		phat, pstd = pp.xing_probs_empirical( Djnk, p = 6, q = 8 )
		row += " & $%0.3f\pm%0.3f$ & $%0.3f\pm%0.3f$ & $%0.3f\pm%0.3f$ & $%0.3f\pm%0.3f$ \\\\ \\hline\\hline\n" % (
		phat[ 0 ], pstd[ 0 ], phat[ 1 ], pstd[ 1 ], phat[ 2 ], pstd[ 2 ], phat[ 3 ], pstd[ 3 ] )
## Commit to table
		body += row 
	return body


body = ""
body += conjecture_table( 'med', 0.5, "" )
body += conjecture_table( 'med', 0.6, "" )
body += conjecture_table( 'med', 0.7, "" )
body += conjecture_table( 'med', 0.8, "" )
body += conjecture_table( 'med', 0.9, "" )
print body


def rejecton_rate( x, alpha = 0.05 ) :
	return 100.0* np.sum( x[ ~np.isnan( x ) ] <= 0.05 ) / np.sum( ~np.isnan( x ), dtype = np.float )

## Generate the body of a table for H = 0.6
def chi_sq_table( method, hurst, body = "" ):
	for kind in kinds :
		row = ""
## get the results
		simd = pp.find_hurst( method, kind, hurst, results )
		if simd is None :
			continue
## Compute teh Chisq_test
		L, H, Njn, Djnk, Vjnde = simd
		stat, dof, pv67 = ss.pooled( Djnk, Njn, p = 6-2, q = 7-2, h = 5 )
		stat, dof, pv78 = ss.pooled( Djnk, Njn, p = 7-2, q = 8-2, h = 5 )
		stat, dof, pv89 = ss.pooled( Djnk, Njn, p = 8-2, q = 9-2, h = 5 )
		stat, dof, pv68 = ss.pooled( Djnk, Njn, p = 6-2, q = 8-2, h = 5 )
		stat, dof, pv79 = ss.pooled( Djnk, Njn, p = 7-2, q = 9-2, h = 5 )
		stat, dof, pv69 = ss.pooled( Djnk, Njn, p = 6-2, q = 9-2, h = 5 )
## Conjectured
		theta = 2**( 1 - 1 / H )
		row += "%s-$%0.2f$" %( kind, H, ) + " & $%0.1f$ & $%0.1f$ & $%0.1f$ & $%0.1f$ & $%0.1f$ & $%0.1f$ \\\\ \\hline \n" % (
		rejecton_rate( pv67 ), rejecton_rate( pv78 ), rejecton_rate( pv89 ),
		rejecton_rate( pv68 ), rejecton_rate( pv79 ), rejecton_rate( pv69 ) )
## Commit to table
		body += row 
	return body

body = ""
body += chi_sq_table( 'med', 0.5, "" )
body += chi_sq_table( 'med', 0.6, "" )
body += chi_sq_table( 'med', 0.7, "" )
body += chi_sq_table( 'med', 0.8, "" )
body += chi_sq_table( 'med', 0.9, "" )
print body






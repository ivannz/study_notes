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

## Returns a list of simulation data
def load_files( files ) :
## The list of loaded results
	results = list( [ ] )
	for file_name in files :
## Parse file name
		prefix, method, run_dttm, size, hurst, replications = re.sub(
			r'(.*)\.npz$', r'\1', os.path.basename( file_name ) ).split( '_' )
## Make the hurst exponent in the filename into a float
		hurst = float( hurst )
## Form a label for the data loaded
		label = prefix + " %.2f" % ( hurst, )
## Load the data
		H, Njn, Djnk, Vjnde = sim_load( file_name, return_durations = False )
		results.append( ( label, hurst, Njn, Djnk, Vjnde ) )
	return results

## Compute the empirical probabilities by averaging across all replications
def xing_probs_empirical( Djnk, p, q ) :
## Fix the levels to analyse
	levels = np.arange( p, q + 1 ) - 1
## Compute the total number of pooled observations
	Cj = np.reshape( np.sum( Djnk[:, levels ], axis = ( 1, 2, ), dtype = np.float ), ( Djnk.shape[ 0 ], 1, ) )
## Compute frequencies of pooled levels
	Pjk = np.where( Cj > 0, np.sum( Djnk[:, levels ], axis = ( 1, ) ) / Cj, 0 )
	return np.average( Pjk, axis = ( 0, ) ), np.std( Pjk, axis = ( 0, ) )

def xing_probs_theoretical( hurst, max_crossings = 20 ) :
## An array of possible sizes of subcrossings
	Z = np.arange( ( max_crossings + 1 ) // 2, dtype = np.int )
## Theoretical distribution
	theta = 2.0 ** ( 1.0 - 1.0 / hurst )
	return 2 * Z + 2, theta * ( 1 - theta ) ** Z

def figure_01( figure, method, kind, p = 3, q = 4 ) :
## Setup view of the polt's main area
	axis = figure.add_subplot( 111 )
## Add a zoomed view of the small size subcrossings region
	#### mini_axis = figure.add_axes( [ .64, .54, .25, .35 ], axisbg = 'w' )
## Setup
	axis.set_xticks( np.arange( 2, 42, 2 ) ) ; axis.grid( )
	#### mini_axis.set_xticks( np.arange( 2, 42, 2 ) ) ; mini_axis.grid( )
	axis.set_color_cycle( plt.cm.rainbow( np.linspace( 0, 1, 5 )[ ::-1 ] ) )
	#### mini_axis.set_color_cycle( plt.cm.rainbow( np.linspace( 0, 1, 5 )[ ::-1 ] ) )
## Draw the empirical probabilities
	for L, H, Njn, Djnk, Vjnde in results[ method ][ kind ] :
## Compute the crossing probabilities
		Phat, Pstd = xing_probs_empirical( Djnk, p = p, q = q )
## Plot theoretical probabilities
		Z, P = xing_probs_theoretical( H, 2 * len( Phat ) )
		axis.plot( Z, P, linestyle = '-', color = 'black' )
		#### mini_axis.plot( Z, P, linestyle = '-', color = 'gray' )
## Plot only positive probabilities
		mask, lwr = Phat > 0, np.where( Phat > Pstd, Pstd, 0 )
## Show 1-sigma error bars
		axis.errorbar( Z[ mask ], Phat[ mask ], yerr = [ lwr[ mask ], Pstd[ mask ] ], fmt = '-o', label = L )
		#### mini_axis.errorbar( Z[ mask ], Phat[ mask ], yerr = [ lwr[ mask ], Pstd[ mask ] ], fmt = '-o', label = L )
## Add proper labels and titles
	axis.set_title( r"""Subcrossing size distribution: empirical against conjectured $\mathbb{P}(Z=2k) = \theta \cdot (1-\theta)^{k-1}$ with $\theta = 2^{1-H^{-1}}$""" )
	axis.set_ylabel( r'Probability (log)' ) ; axis.set_xlabel( r'Subcrossing size' )
## Impose proper geometry
	axis.set_yscale( 'log' ) ; axis.set_ylim( 1e-3, 1.1 ) ; axis.set_xlim( 1.9, 10.1 )
	#### mini_axis.set_yscale( 'log' ) ; mini_axis.set_ylim( .25e-2, 1.1 ) ; mini_axis.set_xlim( 1.9, 6.1 )
## Add legend
	legend = axis.legend( loc = 'lower left', frameon = 1 )
	frame = legend.get_frame( )
	frame.set_facecolor( 'whitesmoke' )
	legend.set_title( r'Processes' )

## For figure 02
def find_hurst( method, kind, hurst ) :
	for L, H, Njn, Djnk, Vjnde in results[ method ][ kind ] :
		if H == hurst :
			return L, H, Njn, Djnk, Vjnde
	return None

def figure_02( figure, method, p, q ) :
## Setup view of the polt's main area
	axis = figure.add_subplot( 111 )
	axis.set_xticks( np.arange( 2, 42, 2 ) ) ; axis.grid( )
## Setup
	colours = plt.cm.rainbow( np.linspace( 0, 1, 5 )[ ::-1 ] )
## Plot H = 0.5
	Z, P = xing_probs_theoretical( 0.5, 2 * 10 )
	axis.plot( Z, P, linestyle = '-', color = 'gray' )
## fBM
	simd = find_hurst( method, 'FBM', 0.5 )
	if simd is not None :
		Phat, Pstd = xing_probs_empirical( simd[ 3 ], p = p, q = q )
		axis.plot( Z[ Phat > 0 ], Phat[ Phat > 0 ], linestyle = '-', color = colours[0], marker = 'o', label = simd[ 0 ] )
## weierstrass
	simd = find_hurst( method, 'WEI', 0.5 )
	if simd is not None :
		Phat, Pstd = xing_probs_empirical( simd[ 3 ], p = p, q = q )
		axis.plot( Z[ Phat > 0 ], Phat[ Phat > 0 ], linestyle = '--', color = colours[0], marker = 'x', label = simd[ 0 ] )
	for hurst, col in zip( [ 0.6, 0.7, 0.8, 0.9, ], colours[ 1: ] ) :
## Plot H = hurst
		Z, P = xing_probs_theoretical( hurst, 2 * 10 )
		axis.plot( Z, P, linestyle = '-', color = 'gray' )
## fBM
		simd = find_hurst( method, 'FBM', hurst )
		if simd is not None :
			Phat, Pstd = xing_probs_empirical( simd[ 3 ], p = p, q = q )
			axis.plot( Z[ Phat > 0 ], Phat[ Phat > 0 ], linestyle = '-', color = col, marker = 'o', label = simd[ 0 ] )
## weierstrass
		simd = find_hurst( method, 'WEI', hurst )
		if simd is not None :
			Phat, Pstd = xing_probs_empirical( simd[ 3 ], p = p, q = q )
			axis.plot( Z[ Phat > 0 ], Phat[ Phat > 0 ], linestyle = '--', color = col, marker = 'x', label = simd[ 0 ] )
## hermite-2
		simd = find_hurst( method, 'HRM-2', hurst )
		if simd is not None :
			Phat, Pstd = xing_probs_empirical( simd[ 3 ], p = p, q = q )
			axis.plot( Z[ Phat > 0 ], Phat[ Phat > 0 ], linestyle = '-.', color = col, marker = '^', label = simd[ 0 ] )
## hermite-3
		simd = find_hurst( method, 'HRM-3', hurst )
		if simd is not None :
			Phat, Pstd = xing_probs_empirical( simd[ 3 ], p = p, q = q )
			axis.plot( Z[ Phat > 0 ], Phat[ Phat > 0 ], linestyle = '-..', color = col, marker = 'v', label = simd[ 0 ] )
## hermite-4
		simd = find_hurst( method, 'HRM-4', hurst )
		if simd is not None :
			Phat, Pstd = xing_probs_empirical( simd[ 3 ], p = p, q = q )
			axis.plot( Z[ Phat > 0 ], Phat[ Phat > 0 ], linestyle = ':', color = col, marker = '>', label = simd[ 0 ] )
## Add proper labels and titles
	axis.set_title( r"""Subcrossing size distribution: empirical against conjectured $\mathbb{P}(Z=2k) = \theta \cdot (1-\theta)^{k-1}$ with $\theta = 2^{1-H^{-1}}$""" )
	axis.set_ylabel( r'Probability (log)' ) ; axis.set_xlabel( r'Subcrossing size' )
## Impose proper geometry
	axis.set_yscale( 'log' ) ; axis.set_ylim( 1e-5, 1.1 ) ; axis.set_xlim( 1.9, 16.1 )
## Add legend
	legend = axis.legend( loc = 2, bbox_to_anchor = ( 0.95, 1 ), borderaxespad = 0.0, frameon = 1, ncol = 1, )
	frame = legend.get_frame( )
	frame.set_facecolor( 'whitesmoke' )
	legend.set_title( r'Processes' )

####################################################################################################
####################################################################################################
####################################################################################################
base_path = os.path.realpath( "./output/final/100-16/" )
process_folders = [ 'FBM_16', 'HRM-2_16-16', 'HRM-3_16-16', 'HRM-4_16-16', 'WEI_16', ]
method_folders = [ 'iqr', 'rng', 'med', ]

results = dict( )
for method in method_folders :
	for process in process_folders :
		path = os.path.join( base_path, process, method )
		if not os.path.exists( path ) :
			continue
## Extract the name of the process from its folder name
		kind = process.split( '_' )[ 0 ]
## Load the results results
		results_method = results.get( method, dict( { } ) )
		results_method[ kind ] = load_files( list_files( path ) )
		results[ method ] = results_method

## With the results loaded, create the first plot
if False :
	kinds = [ 'FBM', 'HRM-2', 'HRM-3', 'HRM-4', 'WEI', ]
	for kind in kinds :
		for method, p, q in [ ('iqr', 6, 7 ), ] : # ('rng', 3, 4), 
			fig = plt.figure( figsize = ( 16, 9 ) )
			figure_01( fig, method, kind, p, q )
			plt.savefig( os.path.join( base_path, 'images', "fig_01_%s_%s.png" % ( method, kind, ) ) )

## FIGURE 02
figure = plt.figure( figsize = ( 16, 9 ) )
figure_02( figure, 'iqr', p = 6, q = 7 )
plt.savefig( os.path.join( base_path, 'images', "fig_02_%s.png" % ( 'iqr', ) ) )


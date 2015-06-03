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
	from operator import itemgetter
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
		H, Njn, Djnk, Vjnde, ( Wjnp, Wbarjn, Wstdjn ) = sim_load( file_name, return_durations = True, without_delta = True )
		results.append( ( label, hurst, Njn, Djnk, Vjnde, Wjnp, Wbarjn, Wstdjn ) )
	return sorted( results, key = itemgetter( 1 ) )

## Compute the empirical probabilities by averaging across all replications
def xing_probs_empirical( Djnk, p, q ) :
## Fix the levels to analyse
	levels = np.arange( p, q + 1 ) - 1
## Compute the total number of pooled observations
	Cj = np.reshape( np.sum( Djnk[:, levels ], axis = ( 1, 2, ), dtype = np.float ), ( Djnk.shape[ 0 ], 1, ) )
## Compute frequencies of pooled levels
	Pjk = np.where( Cj > 0, np.sum( Djnk[:, levels ], axis = ( 1, ) ) / Cj, 0 )
	return np.average( Pjk, axis = ( 0, ) ), np.std( Pjk, axis = ( 0, ) )

## The conjectured crossing probabilities come form a geometric distribution
def xing_probs_theoretical( hurst, max_crossings = 20 ) :
## An array of possible sizes of subcrossings
	Z = np.arange( ( max_crossings + 1 ) // 2, dtype = np.int )
## Theoretical distribution
	theta = 2.0 ** ( 1.0 - 1.0 / hurst )
	return 2 * Z + 2, theta * ( 1 - theta ) ** Z

## Return the empirical distribution of excursions
def xcur_probs_empirical( Vjnde, p, q ) :
	levels = np.arange( p, 1 + q ) - 1
## Compute the pooled number of excursions of both types conditional on the crossing direction.
	Cjd = np.reshape( np.sum( Vjnde[:, levels ], axis = ( 1, 3, ), dtype = np.float ),
		( ( Vjnde.shape[ 0 ], Vjnde.shape[ 2 ], 1 ) ) )
## Estimate the chances of excursions conditional on the crossing direction.
	Pjde = np.where( Cjd > 0, np.sum( Vjnde[:, levels ], axis = ( 1, ) ) / Cjd, 0 )
## We are interested in the diagonal entries: +-|++ and -+|--.
	return np.diagonal( Pjde, offset = 0, axis1 = 1, axis2 = 2 )

## Estimate the hurst exponent from the mean number of subcrossings 
def xing_hurst( Djnk, p, q ) :
	levels = np.arange( p, 1 + q ) - 1
## Get the number of 
	Cjn = np.reshape( np.sum( Djnk[ :, levels ], axis = ( 2, ), dtype = np.float ), ( Djnk.shape[ 0 ], q - p + 1, 1, ) )
	Pjnk = np.where( Cjn > 0, Djnk[ :, levels ] / Cjn, 0 )
	Mjn = 2 * np.tensordot( Pjnk, 1 + np.arange( Pjnk.shape[ -1 ] ), axes = ( 2, 0 ) )
	Hjn = np.log( 2 ) / np.log( Mjn )
	return np.average( Hjn, axis = ( 0, ) ).flatten( ), np.std( Hjn, axis = ( 0, ) ).flatten( )

def xing_geom_mle( Djnk, p, q ) :
## Fix the levels to analyse
	levels = np.arange( p, q + 1 ) - 1
## Compute the total number of pooled observations
	Cj = np.reshape( np.sum( Djnk[ :, levels, :-1 ], axis = ( 1, 2, ), dtype = np.float ), Djnk.shape[ :1 ] + ( 1, ) )
## Compute frequencies of pooled levels
	Pjk = np.where( Cj > 0, np.sum( Djnk[:, levels ], axis = ( 1, ) ) / Cj, 0 )
## Compute the truncated mean
	Mj = np.tensordot( Pjk, np.arange( Pjk.shape[ -1 ], dtype = np.float ), axes = ( 1, 0 ) )
	return 1.0 / ( 1 + Mj )

def xing_atomic_geom_mle( Djnk, p, q ) :
## Fix the levels to analyse
	levels = np.arange( p, q + 1 ) - 1
## Compute the total number of pooled observations
	Cj = np.reshape( np.sum( Djnk[ :, levels, 1:-1 ], axis = ( 1, 2, ), dtype = np.float ), Djnk.shape[ :1 ] + ( 1, ) )
## Compute frequencies of pooled levels
	Pjk = np.where( Cj > 0, np.sum( Djnk[:, levels, 1: ], axis = ( 1, ) ) / Cj, 0 )
## Compute the truncated mean
	Mj = np.tensordot( Pjk, np.arange( Pjk.shape[ -1 ], dtype = np.float ), axes = ( 1, 0 ) )
## Get theta
	theta_j = 1.0 / ( 1 + Mj )
## Compute delta
	N2j = np.sum( Djnk[ :, levels, :1 ], axis = ( 1, 2, ), dtype = np.float )
	Cj = np.sum( Djnk[ :, levels ], axis = ( 1, 2, ), dtype = np.float )
	delta_j = np.where( Cj > 0, N2j / Cj - theta_j, 0 )
## Get delta
	return theta_j, delta_j

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
	for L, H, Njn, Djnk, Vjnde, Wjnp, Wbarjn, Wstdjn in results[ method ][ kind ] :
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
def find_hurst( method, kind, hurst, results ) :
	if method not in results:
		return None
	if kind not in results[ method ] :
		return None
	for L, H, Njn, Djnk, Vjnde, Wjnp, Wbarjn, Wstdjn in results[ method ][ kind ] :
		if H == hurst :
			return L, H, Njn, Djnk, Vjnde, Wjnp, Wbarjn, Wstdjn
	return None

def figure_02( figure, method, p, q ) :
## Setup view of the polt's main area
	axis = figure.add_subplot( 111 )
	axis.set_xticks( np.arange( 2, 42, 2 ) ) ; axis.grid( )
## Setup
	colours = plt.cm.rainbow( np.linspace( 0, 1, 5 )[ ::-1 ] )
## Plot H = 0.5
	Z, P = xing_probs_theoretical( 0.5, 2 * 20 )
	axis.plot( Z, P, linestyle = '-', color = 'gray' )
## fBM
	simd = find_hurst( method, 'FBM', 0.5, results )
	if simd is not None :
		Phat, Pstd = xing_probs_empirical( simd[ 3 ], p = p, q = q )
		axis.plot( Z[ np.where( Phat > 0 )[ 0 ] ], Phat[ Phat > 0 ], linestyle = '-', color = colours[0], marker = 'o', label = simd[ 0 ] )
## weierstrass
	simd = find_hurst( method, 'WEI', 0.5, results )
	if simd is not None :
		Phat, Pstd = xing_probs_empirical( simd[ 3 ], p = p, q = q )
		axis.plot( Z[ np.where( Phat > 0 )[ 0 ] ], Phat[ Phat > 0 ], linestyle = ':', color = colours[0], marker = '*', markersize = 10, label = simd[ 0 ] )
	for hurst, col in zip( [ 0.6, 0.7, 0.8, 0.9, ], colours[ 1: ] ) :
## Plot H = hurst
		Z, P = xing_probs_theoretical( hurst, 2 * 20 )
		axis.plot( Z, P, linestyle = '-', color = 'gray' )
## fBM
		simd = find_hurst( method, 'FBM', hurst, results )
		if simd is not None :
			Phat, Pstd = xing_probs_empirical( simd[ 3 ], p = p, q = q )
			axis.plot( Z[ np.where( Phat > 0 )[ 0 ] ], Phat[ Phat > 0 ], linestyle = '-', color = col, marker = 'o', label = simd[ 0 ] )
## weierstrass
		simd = find_hurst( method, 'WEI', hurst, results )
		if simd is not None :
			Phat, Pstd = xing_probs_empirical( simd[ 3 ], p = p, q = q )
			axis.plot( Z[ np.where( Phat > 0 )[ 0 ] ], Phat[ Phat > 0 ], linestyle = ':', color = col, marker = '*', markersize = 10, label = simd[ 0 ] )
## hermite-2
		simd = find_hurst( method, 'HRM-2', hurst, results )
		if simd is not None :
			Phat, Pstd = xing_probs_empirical( simd[ 3 ], p = p, q = q )
			axis.plot( Z[ np.where( Phat > 0 )[ 0 ] ], Phat[ Phat > 0 ], linestyle = '--', color = col, marker = '^', label = simd[ 0 ] )
## hermite-3
		simd = find_hurst( method, 'HRM-3', hurst, results )
		if simd is not None :
			Phat, Pstd = xing_probs_empirical( simd[ 3 ], p = p, q = q )
			axis.plot( Z[ np.where( Phat > 0 )[ 0 ] ], Phat[ Phat > 0 ], linestyle = '--', color = col, marker = 'v', label = simd[ 0 ] )
## hermite-4
		simd = find_hurst( method, 'HRM-4', hurst, results )
		if simd is not None :
			Phat, Pstd = xing_probs_empirical( simd[ 3 ], p = p, q = q )
			axis.plot( Z[ np.where( Phat > 0 )[ 0 ] ], Phat[ Phat > 0 ], linestyle = '--', color = col, marker = 's', label = simd[ 0 ] )
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

## The main plotter for the third figure
def __plot_excursions( axis, method, hurst, p, q, direction = 0 ) :
	axis.set_ylim( 0.35, 0.85 )
	# axis.set_ylim( -0.1, 1.1 )
	axis.spines[ 'top' ].set_visible( False ) ; axis.spines[ 'bottom' ].set_visible( False )
	axis.xaxis.set_ticks_position( 'none' ) ; axis.yaxis.set_ticks_position( 'left' )
	axis.axhline( y = 1.0 / np.sqrt( 2 ** ( 1.0 / hurst ) ), linewidth = 2, color = 'black' )
	boxes, labels = list( ), list( )
	for kind in [ 'FBM', 'WEI', 'HRM-2', 'HRM-3', 'HRM-4', ] :
		simd = find_hurst( method, kind, hurst, results )
		if simd is not None :
			Pje = xcur_probs_empirical( simd[ 4 ], p, q )
			boxes.append( Pje[:, direction ] )
			labels.append( simd[ 0 ] )
	if boxes :
		axis.boxplot( boxes )
		axis.set_xticklabels( labels, rotation = -90, position = (0,0.1) )

def figure_03( figure, method, p, q, direction ) :
	for hurst, subplot in [  ( 0.5, 1 ), ( 0.6, 2 ), ( 0.7, 3 ), ( 0.8, 4 ), ( 0.9, 5 ), ] :
		__plot_excursions( figure.add_subplot( 1, 5, subplot ), method, hurst, p = p, q = q, direction = direction )

def figure_04( figure, method, kind ) :
	axis = figure.add_subplot( 111 )
	# colours = plt.cm.rainbow( np.linspace( 0, 1, 5 )[ ::-1 ] )
	axis.set_color_cycle( plt.cm.rainbow( np.linspace( 0, 1, 5 )[::-1] ) )
## Get the highest level plotted
	highest_level = -1
	for L, H, Njn, Djnk, Vjnde, Wjnp, Wbarjn, Wstdjn in results[ method ][ kind ] :
		Hhat, Hstd = xing_hurst( Djnk, 1, Djnk.shape[ 1 ] )
		axis.axhline( y = H, color = 'black', linestyle = '-.' )
		lwr, mask = np.where( Hhat > Hstd, Hstd, 0 ), Hhat > H * 0.75
		Ln = 1 + np.arange( len( Hhat ) )
		axis.errorbar( Ln[ mask ], Hhat[ mask ], yerr = [ lwr[ mask ], Hstd[ mask ] ],
			fmt = '-o', label = L )
		highest_level = max( np.max( Ln[ mask ] ), highest_level )
	axis.set_xlim( 0.9, highest_level + 1 ) ; axis.set_ylim( 0.45, 1.01 )
	## Add a legend with white opaque background.
	axis.set_title( 'Crossing tree estimates of the Hurst exponent' )
	axis.set_xlabel( 'Level of the crossing tree' ) ; axis.set_ylabel( r"""$\hat{H}$""", rotation = 0 )
	legend = axis.legend( loc = 'lower right', frameon = 1 )
	frame = legend.get_frame( )
	frame.set_facecolor( 'whitesmoke' )

def figure_05( figure, method, p, q ) :
## Setup view of the polt's main area
	axis = figure.add_subplot( 111 )
## Setup
	colours = plt.cm.rainbow( np.linspace( 0, 1, 5 )[ ::-1 ] )
## Fix the levels
	levels = np.arange( p, q + 1 )
	highest_level = -1
## Plot H = 0.5
	axis.axhline( y = 0.5, color = 'gray', linestyle = '-' )
## fBM
	simd = find_hurst( method, 'FBM', 0.5, results )
	if simd is not None :
		Hhat, Hstd = xing_hurst( simd[ 3 ], p, q )
		mask = np.where( Hhat > 0.80 * 0.5 )[ 0 ]
		highest_level = max( np.max( levels[ mask ] ), highest_level )
		axis.plot( levels[ mask ], Hhat[ mask ], linestyle = '-', color = colours[ 0 ], marker = 'o', label = simd[ 0 ] )
## weierstrass
	simd = find_hurst( method, 'WEI', 0.5, results )
	if simd is not None :
		Hhat, Hstd = xing_hurst( simd[ 3 ], p, q )
		mask = np.where( Hhat > 0.80 * 0.5 )[ 0 ]
		highest_level = max( np.max( levels[ mask ] ), highest_level )
		axis.plot( levels[ mask ], Hhat[ mask ], linestyle = ':', color = colours[ 0 ], marker = '*', markersize = 10, label = simd[ 0 ] )
	for hurst, col in zip( [ 0.6, 0.7, 0.8, 0.9, ], colours[ 1: ] ) :
## Plot H = hurst
		axis.axhline( y = hurst, color = 'gray', linestyle = '-' )
## fBM
		simd = find_hurst( method, 'FBM', hurst, results )
		if simd is not None :
			Hhat, Hstd = xing_hurst( simd[ 3 ], p, q )
			mask = np.where( Hhat > 0.80 * hurst )[ 0 ]
			highest_level = max( np.max( levels[ mask ] ), highest_level )
			axis.plot( levels[ mask ], Hhat[ mask ], linestyle = '-', color = col, marker = 'o', label = simd[ 0 ] )
## weierstrass
		simd = find_hurst( method, 'WEI', hurst, results )
		if simd is not None :
			Hhat, Hstd = xing_hurst( simd[ 3 ], p, q )
			mask = np.where( Hhat > 0.80 * hurst )[ 0 ]
			highest_level = max( np.max( levels[ mask ] ), highest_level )
			axis.plot( levels[ mask ], Hhat[ mask ], linestyle = ':', color = col, marker = '*', markersize = 10, label = simd[ 0 ] )
## hermite-2
		simd = find_hurst( method, 'HRM-2', hurst, results )
		if simd is not None :
			Hhat, Hstd = xing_hurst( simd[ 3 ], p, q )
			mask = np.where( Hhat > 0.80 * hurst )[ 0 ]
			highest_level = max( np.max( levels[ mask ] ), highest_level )
			axis.plot( levels[ mask ], Hhat[ mask ], linestyle = '--', color = col, marker = '^', label = simd[ 0 ] )
## hermite-3
		simd = find_hurst( method, 'HRM-3', hurst, results )
		if simd is not None :
			Hhat, Hstd = xing_hurst( simd[ 3 ], p, q )
			mask = np.where( Hhat > 0.80 * hurst )[ 0 ]
			highest_level = max( np.max( levels[ mask ] ), highest_level )
			axis.plot( levels[ mask ], Hhat[ mask ], linestyle = '--', color = col, marker = 'v', label = simd[ 0 ] )
## hermite-4
		simd = find_hurst( method, 'HRM-4', hurst, results )
		if simd is not None :
			Hhat, Hstd = xing_hurst( simd[ 3 ], p, q )
			mask = np.where( Hhat > 0.80 * hurst )[ 0 ]
			highest_level = max( np.max( levels[ mask ] ), highest_level )
			axis.plot( levels[ mask ], Hhat[ mask ], linestyle = '--', color = col, marker = 's', label = simd[ 0 ] )
	axis.set_xticks( np.arange( p, highest_level, 1 ) ) ; axis.grid( )
	axis.set_xlim( p - 0.1, highest_level + 1.1 ) ; axis.set_ylim( 0.45, 1.01 )
## Add proper labels and titles
	axis.set_title( r"""Crossing tree estimates of the Hurst exponent: $H$-SSSI processes""" )
	axis.set_ylabel( r"""$\hat{H}$""" ) ; axis.set_xlabel( r'Level of the crossing tree' )
## Add legend
	legend = axis.legend( loc = 2, bbox_to_anchor = ( 0.95, 1 ), borderaxespad = 0.0, frameon = 1, ncol = 1, )
	frame = legend.get_frame( )
	frame.set_facecolor( 'whitesmoke' )
	legend.set_title( r'Processes' )

## The main plotter for the third figure
def __plot_mle_probs( axis, method, hurst, p, q, with_atom = False ) :
	axis.set_ylim( -0.1, 1.1 )
	axis.spines[ 'top' ].set_visible( False ) ; axis.spines[ 'bottom' ].set_visible( False )
	axis.xaxis.set_ticks_position( 'none' ) ; axis.yaxis.set_ticks_position( 'left' )
	axis.axhline( y = 2.0 ** ( 1.0 - 1.0 / hurst ), linewidth = 2, color = 'black' )
	boxes, labels = list( ), list( )
	for kind in [ 'FBM', 'WEI', 'HRM-2', 'HRM-3', 'HRM-4', ] :
		simd = find_hurst( method, kind, hurst, results )
		if simd is not None :
			if not with_atom :
				theta_j = xing_geom_mle( simd[ 3 ], p, q )
			else :
				theta_j, delta_j = xing_atomic_geom_mle( simd[ 3 ], p, q )
			boxes.append( theta_j )
			labels.append( simd[ 0 ] )
	if boxes :
		axis.boxplot( boxes )
		axis.set_xticklabels( labels, rotation = -90, position = (0,0.2) )

def figure_06( figure, method, p, q ) :
	for hurst, subplot in [  ( 0.5, 1 ), ( 0.6, 2 ), ( 0.7, 3 ), ( 0.8, 4 ), ( 0.9, 5 ), ] :
		__plot_mle_probs( figure.add_subplot( 1, 5, subplot ), method, hurst, p = p, q = q, with_atom = False )

def figure_07( figure, method, p, q ) :
	for hurst, subplot in [  ( 0.5, 1 ), ( 0.6, 2 ), ( 0.7, 3 ), ( 0.8, 4 ), ( 0.9, 5 ), ] :
		__plot_mle_probs( figure.add_subplot( 1, 5, subplot ), method, hurst, p = p, q = q, with_atom = True )


## Show mean crossing durations and their standard error
def figure_08( figure, method, kind ) :
	figure = plt.figure( figsize = ( 16, 9 ) )
	axis = figure.add_subplot( 111 )
	axis.set_color_cycle( plt.cm.rainbow( np.linspace( 0, 1, 5 )[ ::-1 ] ) )
	for L, H, Njn, Djnk, Vjnde, Wjnp, Wbarjn, Wstdjn in results[ method ][ kind ] :
		Z = np.arange( Wbarjn.shape[ 1 ] ) + 1
		scaled_wbarjn = Wbarjn * ( 2.0 ** ( - np.arange( Wbarjn.shape[ 1 ] ).reshape( ( 1, -1, 1 ) ) / H ) )
		What, Wstd = np.average( scaled_wbarjn, axis = ( 0, ) ), np.std( scaled_wbarjn, axis = ( 0, ) )
## Plot only positive probabilities
		mask, lwr = What[:,0] > 0, np.where( What > Wstd, Wstd, 0 )
## Show 1-sigma error bars
		axis.errorbar( Z[ mask ], What[ mask ], yerr = [ lwr[ mask ], Wstd[ mask ] ], fmt = '-', label = L )
		# axis.boxplot( [ Wbarjn[:,n] * 2 ** ( - n / H ) for n in range( Wbarjn.shape[ 1 ] ) ] )
	axis.set_xlabel( 'Level of the crossing tree' )
	axis.set_ylabel( r"""$2^{-n H^{-1}}\mathbb{E}{W}_n$""", rotation = 90 )
## Impose proper geometry
	axis.set_yscale( 'log' ) ; axis.set_xlim( 0.9, 20.1 )
## Add legend
	legend = axis.legend( loc = 'lower left', frameon = 1 )
	frame = legend.get_frame( )
	frame.set_facecolor( 'whitesmoke' )
	legend.set_title( r'Processes' )

## For figure 02

####################################################################################################
####################################################################################################
####################################################################################################
if __name__ == '__main__' :
	base_path = os.path.realpath( "./output/10000-17" )
	process_folders = [ 'FBM_17', 'HRM-2_17-16', 'HRM-3_17-16', 'HRM-4_17-16', 'WEI_17', ]
	# base_path = os.path.realpath( "./output/1000-21" )
	# process_folders = [ 'FBM_21', ]
	# base_path = os.path.realpath( "./output/release/1000-16" )
	# process_folders = [ 'FBM_16', 'HRM-2_16-16', 'HRM-3_16-16', 'HRM-4_16-16', 'WEI_16', ]


	## Quickly extract the info
	kinds = [ process.split( '_' )[ 0 ] for process in process_folders ]

	## Method name and the tree levels to use (chosen heuristically).
	# methods = [ ( 'med', 7, 8 ), ( 'iqr', 6, 7 ), ('rng', 3, 4 ), ]
	methods = [ ( 'med', 7, 8 ), ]

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
			results_method[ kind ] = load_files( list_files( path ) )
			results[ method ] = results_method

	## FIGURE 01
	if False :
		for method, p, q in methods :
			if method not in results :
				continue
			for kind in kinds :
				if kind not in results[ method ] :
					continue
				fig = plt.figure( figsize = ( 16, 9 ) )
				figure_01( fig, method, kind, p, q )
				plt.savefig( os.path.join( base_path, 'pdf', "fig_01_%s_%s.pdf" % ( method, kind, ) ) , format = 'pdf' )

	## FIGURE 02
	if True :
		for method, p, q in methods :
			if method in results :
				figure = plt.figure( figsize = ( 16, 9 ) )
				figure_02( figure, method, p = p, q = q )
				plt.savefig( os.path.join( base_path, 'pdf', "fig_02_%s.pdf" % ( method, ) ) , format = 'pdf' )

	## FIGURE 03
	if False :
		for direction, name in [ (0, 'up-down'), (1, 'down-up') ] :
			for method, p, q in methods :
				if method in results :
					figure = plt.figure( figsize = ( 16, 8 ) )
					plt.subplots_adjust( hspace = 0.4 )
					if direction == 0 :
						figure.suptitle( r"Empirical probability of an %s excursion conditional on an up-crossing" % ( name, ) )
					else :
						figure.suptitle( r"Empirical probability of an %s excursion conditional on a down-crossing" % ( name, ) )
					figure_03( figure, method, p = p, q = q, direction = direction )
					plt.savefig( os.path.join( base_path, 'pdf', "fig_03_%s_%s.pdf" % ( name, method, ) ) , format = 'pdf' )

	## FIGURE 04
	if False :
		for method, p, q in methods :
			if method not in results :
				continue
			for kind in kinds :
				if kind not in results[ method ] :
					continue
				figure = plt.figure( figsize = ( 16, 9 ) )
				figure_04( figure, method, kind )
				plt.savefig( os.path.join( base_path, 'pdf', "fig_04_%s_%s.pdf" % ( method, kind, ) ) , format = 'pdf' )

	## FIGURE 05
	if False :
		for method, p, q in methods :
			if method not in results :
				continue
			figure = plt.figure( figsize = ( 16, 9 ) )
			figure_05( figure, method, 1, 20 )
			plt.savefig( os.path.join( base_path, 'pdf', "fig_05_%s.pdf" % ( method, ) ) , format = 'pdf' )

	## FIGURE 06
	if False :
		for method, p, q in methods :
			if method in results :
				figure = plt.figure( figsize = ( 16, 5 ) )
				plt.subplots_adjust( hspace = 0.4 )
				figure.suptitle( r"""MLE estimate of the $Z_k \sim$Geom$(\theta)$""" )
				figure_06( figure, method, p = p, q = q )
				plt.savefig( os.path.join( base_path, 'pdf', "fig_06_%s.pdf" % ( method, ) ) , format = 'pdf' )

	## FIGURE 07
	if False :
		for method, p, q in methods :
			if method in results :
				figure = plt.figure( figsize = ( 16, 5 ) )
				plt.subplots_adjust( hspace = 0.4 )
				figure.suptitle( r"""MLE estimate of the $Z_k \sim$Geom$(\theta)$ with atom at $2$""" )
				figure_07( figure, method, p = p, q = q )
				plt.savefig( os.path.join( base_path, 'pdf', "fig_07_%s.pdf" % ( method, ) ) , format = 'pdf' )

	## FIGURE 08
	if False :
		for method, p, q in methods :
			if method not in results :
				continue
			for kind in kinds :
				if kind not in results[ method ] :
					continue
				fig = plt.figure( figsize = ( 16, 9 ) )
				figure.suptitle( r"""Average scaled crossing durations for fBm""" )
				figure_08( fig, method, kind )
				plt.savefig( os.path.join( base_path, 'pdf', "fig_08_%s_%s.pdf" % ( method, kind, ) ) , format = 'pdf' )

	if False :
		Phat, Pstd = xing_probs_empirical( simd[ 4 ], p = p, q = q )
		axis.plot( Z[ np.where( Phat > 0 )[ 0 ] ], Phat[ Phat > 0 ], linestyle = '--', color = col, marker = 's', label = simd[ 0 ] )
		figure = plt.figure( figsize = ( 16, 9 ) )
		axis = figure.add_subplot( 111 )
		boxes, labels = list( ), list( )
		for L, H, Njn, Djnk, Vjnde, Wjnp, Wbarjn, Wstdjn in results[ 'med' ][ 'WEI' ] :
			Pje = xcur_probs_empirical( Vjnde, p, q )
			boxes.append( Pje[:,0] )
			labels.append( L )
		axis.boxplot( boxes, labels = labels )
		plt.show( figure )


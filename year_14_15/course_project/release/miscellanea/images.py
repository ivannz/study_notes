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

def __get_xing_probs( Djnk, H, p = 5, q = 7 ) :
## Fix the levels to analyse
	levels = np.arange( p, q + 1 ) - 1
## Compute the total number of pooled observations
	Cj = np.reshape( np.sum( Djnk[:, levels ], axis = ( 1, 2, ), dtype = np.float ), ( Djnk.shape[ 0 ], 1, ) )
## Compute frequencies of pooled levels
	Pjk = np.where( Cj > 0, np.sum( Djnk[:, levels ], axis = ( 1, ) ) / Cj, 0 )
## An array of possible sizes of subcrossings
	Zk = np.reshape( 1 + np.arange( Pjk.shape[ -1 ] ), ( 1, Pjk.shape[ -1 ] ) )
## Theoretical distribution
	theta = 2.0**( 1.0 - 1.0 / H )
	Pk = theta * ( 1 - theta )**( Zk - 1 )
## Return
	return 2 * Zk.flatten(), Pk.flatten(), np.average( Pjk, axis = ( 0, ) ), np.std( Pjk, axis = ( 0, ) )

def __draw_probs( ax, sim_data, p = 2, q = 3, theor_color = 'gray' ) :
	ax.set_xticks( np.arange( 2, 42, 2 ) ) ; ax.grid( )
	ax.set_color_cycle( plt.cm.rainbow( np.linspace( 0, 1, max( len( files ), 10 ) )[ ::-1 ] ) )
	for L, H, Njn, Djnk, Vjnde in sim_data :
## Compute the crossing probabilities
		Zk, Pk, Phat, Pstd = __get_xing_probs( Djnk, H, p = p, q = q )
## PLot reference values
		ax.plot( Zk, Pk, linestyle = '-', color = theor_color )
## Plot only positive probabilities
		mask, lwr = Phat > 0, np.where( Phat > Pstd, Pstd, 0 )
## Show 1-sigma error bars
		ax.errorbar( Zk[ mask ], Phat[ mask ], yerr = [ lwr[ mask ], Pstd[ mask ] ], fmt = '-o', label = L )

def draw_probs( sim_data ) :
## Create an appropriately sized figure
	fig = plt.figure( figsize = ( 16, 9 ) )
	ax = plt.subplot( 111 )
## Setup view of the polt's main area
	__draw_probs( ax, sim_data, theor_color = 'black' )
	ax.set_yscale( 'log' ) ; ax.set_ylim( 1e-3, 1 ) ; ax.set_xlim( 2, 10 )
## Add proper labels and titles
	ax.set_title( r"""Subcrossing size distribution: empirical against conjectured $\mathbb{P}(Z=2k) = \theta \cdot (1-\theta)^{k-1}$ with $\theta = 2^{1-H^{-1}}$""" )
	ax.set_ylabel( r'Probability (log)' ) ; ax.set_xlabel( r'Subcrossing size' )
## Add a zoomed view of the small size subcrossings region
	ax2 = plt.axes( [ .64, .54, .25, .35 ], axisbg = 'w' )
	__draw_probs( ax2, sim_data, theor_color = 'gray' )
	ax2.set_yscale( 'log' ) ; ax2.set_ylim( .25e-2, 1 ) ; ax2.set_xlim( 1, 7 )
## Add legend
	legend = ax.legend( loc = 'lower left', frameon = 1 )
	frame = legend.get_frame( )
	frame.set_facecolor( 'whitesmoke' )
	legend.set_title( r'Empirical' )
## Output
	return fig

## Now do the hurst plots
def xing_hurst( file_name, p = 5, q = 7 ) :
	levels = np.arange( p, q + 1 ) - 1
	H, Njn, Djnk, Vjnde = sim_load( file_name )
	Cjn = np.reshape( np.sum( Djnk, axis = ( 2, ), dtype = np.float ), Djnk.shape[ :-1 ] + ( 1, ) )
	Zk = 2 * np.reshape( 1 + np.arange( Djnk.shape[ -1 ] ), ( 1, 1, Djnk.shape[ -1 ] ) )
	Pjnk = np.where( Cjn > 0, Djnk / Cjn, 0 )
	Mjn = np.reshape( np.sum( Pjnk * Zk, axis = ( 2, ) ), Cjn.shape )
	Hjn = np.log( 2 ) / np.log( Mjn )
	return ( ( 1 + np.arange( Djnk.shape[ -2 ] ) )[levels], H,
		np.average( Hjn, axis = (0, ) ).flatten( )[levels],
		np.std( Hjn, axis = (0, ) ).flatten( )[levels] )


# files = list_files( './output/recent/HRM-3_20-16/iqr/' )
files = list_files( './output/releae/HRM-2_16-16/rng/' )
results = load_files( files )
fig_01 = draw_probs( results )
plt.savefig( './output/plots/HRM-2_16-16-rng-log_probs.png', fig )







p, q = 6, 8
method = 'med'
kind = 'FBM'

L, H, Njn, Djnk, Vjnde = results[ method ][ kind ][ 4 ]
tj, dj = xing_atomic_geom_mle( Djnk, p, q )
gj = xing_geom_mle( Djnk, p, q )
axis = plt.subplot( 111 )
axis.axhline( y = 2.0 ** ( 1.0 - 1.0 / H ), linewidth = 2, color = 'black' )
axis.boxplot( [ gj, tj, dj ] )
plt.show( )


plt.figure( figsize = ( 16, 9 ) )
plt.subplot( 131 )
stat, dof, pv = ss.pooled( Djnk, Njn, p = 6, q = 7, h = 5 )
plt.hist( pv[ ~ np.isnan( pv ) ] )
plt.subplot( 132 )
stat, dof, pv = ss.pooled( Djnk, Njn, p = 7, q = 8, h = 5 )
plt.hist( pv[ ~ np.isnan( pv ) ] )
plt.subplot( 133 )
stat, dof, pv = ss.pooled( Djnk, Njn, p = 6, q = 8, h = 5 )
plt.hist( pv[ ~ np.isnan( pv ) ] )
plt.show( )

plt.figure( figsize = ( 16, 9 ) )
plt.subplot( 331 )
stat, dof, pv = ss.pooled( Djnk, Njn, p = 6, q = 5+2, h = 6+2 )
plt.hist( pv[ ~ np.isnan( pv ) ] )
plt.subplot( 332 )
stat, dof, pv = ss.pooled( Djnk, Njn, p = 6, q = 6+2, h = 6+2 )
plt.hist( pv[ ~ np.isnan( pv ) ] )
plt.subplot( 333 )
stat, dof, pv = ss.pooled( Djnk, Njn, p = 6, q = 7+2, h = 6+2 )
plt.hist( pv[ ~ np.isnan( pv ) ] )
plt.subplot( 334 )
stat, dof, pv = ss.pooled( Djnk, Njn, p = 6, q = 5+2, h = 5+2 )
plt.hist( pv[ ~ np.isnan( pv ) ] )
plt.subplot( 335 )
stat, dof, pv = ss.pooled( Djnk, Njn, p = 6, q = 6+2, h = 5+2 )
plt.hist( pv[ ~ np.isnan( pv ) ] )
plt.subplot( 336 )
stat, dof, pv = ss.pooled( Djnk, Njn, p = 6, q = 7+2, h = 5+2 )
plt.hist( pv[ ~ np.isnan( pv ) ] )
plt.subplot( 337 )
stat, dof, pv = ss.pooled( Djnk, Njn, p = 6, q = 5+2, h = 4+2 )
plt.hist( pv[ ~ np.isnan( pv ) ] )
plt.subplot( 338 )
stat, dof, pv = ss.pooled( Djnk, Njn, p = 6, q = 6+2, h = 4+2 )
plt.hist( pv[ ~ np.isnan( pv ) ] )
plt.subplot( 339 )
stat, dof, pv = ss.pooled( Djnk, Njn, p = 6, q = 7+2, h = 4+2 )
plt.hist( pv[ ~ np.isnan( pv ) ] )
plt.show( )

plt.figure( figsize = ( 16, 9 ) )
plt.subplot( 221 )
stat, dof, pv = ss.pooled( Djnk, Njn, p = 6, q = 7, h = 7 )
plt.hist( pv[ ~ np.isnan( pv ) ] )
plt.subplot( 222 )
stat, dof, pv = ss.pooled( Djnk, Njn, p = 5, q = 6, h = 7 )
plt.hist( pv[ ~ np.isnan( pv ) ] )
plt.subplot( 223 )
stat, dof, pv = ss.pooled( Djnk, Njn, p = 5, q = 8, h = 7 )
plt.hist( pv[ ~ np.isnan( pv ) ] )
plt.subplot( 224 )
stat, dof, pv = ss.pooled( Djnk, Njn, p = 6, q = 8, h = 7 )
plt.hist( pv[ ~ np.isnan( pv ) ] )
plt.show( )


plt.figure( figsize = ( 16, 9 ) )
plt.subplot( 121 )
stat, dof, pv = ss.pooled( Djnk, Njn, p = 6, q = 7, h = 7 )
plt.hist( pv[ ~ np.isnan( pv ) ] )
plt.subplot( 122 )
stat, dof, pv = ss.pooled( Djnk, Njn, p = 5, q = 6, h = 7 )
plt.hist( pv[ ~ np.isnan( pv ) ] )
plt.show( )




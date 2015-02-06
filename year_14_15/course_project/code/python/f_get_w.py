# coding: utf-8
import numpy as np

# times, points and subcrossings for continuous process
#
# x: signal to be analysed
# t: time ( mandatory )
# levels: for each level in levels crossings will be calculated at scale delta * 2^level
# delta: the base scale
# deleteFirst: if 1 then delete first crossing at each level
# 
# w: crossing duration ( cell structure, size of the number of levels )
# subx : number of subcrossings ( cell structure, size of the number of levels, first empty )
# hit_point: cell struct., hitting levels 
# hit_time: cell struct., hitting times
# Y.Shen 12.2002, O.D.Jones 7.2003, 6.2007
# P - O. Amblard & G. Decrouez 2011
# Nazarov 2015
#
#  [w, subx, hit_point, hit_time] = f_get_w( fname, levels, delta, deleteFirst )
## X, T -- aligned value and time series of the process
##  X, T = ( {x(t_i)}_{i=1}^n, {t_i}_{i=1}^n )
def f_get_w( T, X, levels, delta, deleteFirst = False ) :
## Rescale the sample path, so that the grid base scale is 1.0 
	ht, hp = grid_hitting( T, ( X - X[ 0 ] ) / delta, levels, deleteFirst, X[ 0 ], delta )
	subxing = get_subxings( ht, hp )
## Get the waiting times and rescale the hitting points
	wt = [ np.diff( lht ) for lht in ht ]
	return ( ht, hp, subxing, wt )

## For a given set of grid levels compute hitting times
def grid_hitting( T, X, levels, skipFirst = False, mu = 0.0, sigma = 1.0 ) :
	hp = list() ; ht = list()
	for level in levels :
		lht, lhp = grid_hits( T, X, scale = 2 ** level, skipFirst = skipFirst, mu = mu, sigma = sigma )
		ht.append( lht )
		hp.append( lhp )
	return ( ht, hp )

## Compute the crossings by a process X_T of the uniformly
##  spaced grid.
def grid_hits( T, X, scale, skipFirst = False, mu = 0.0, sigma = 1.0 ) :
	lht = list([ 0 ]) ; lhp = list([ mu ])
## Assume the process starts at zero, so intitally crossings
##  of the the zero-th grid level can be easily ignored.
	last_hit = 0
	for t in xrange( 1, len( X ) ) :
		if X[ t - 1 ] == X[ t ] : continue
## Check how the process has moved.
		if X[ t - 1 ] < X[ t ] :
## The process moved in upward direction:
			direction  = 1
## Find the first grid level it should have crossed
			level0 = np.ceil(  X[ t - 1 ] / scale )
## Find the last level, it should have passed on its way up 
			level1 = np.floor( X[ t     ] / scale )
## If process hasn't made a crossing, that is the endpoints
##  are between two levels but within one band, then 
		else :
			direction  = -1
			level0 = np.floor( X[ t - 1 ] / scale )
			level1 = np.ceil(  X[ t     ] / scale )
## For each grid level between the first and the last levels
##  crossed by the process
##		print level0, level1, direction
		for j in np.arange( level0, level1 + direction, direction ) :
## It is not a crossing if the process returned to the level
##  of the previous crossing
##  T^n_{k+1} \defn \inf\{  . t > T^n_k : \, X_t\in \partial B_n,\, X_t \neq X_{T^n_k} \}
			if j == last_hit : continue
			lhp.append( j * scale * sigma + mu )
## Due to discrete data the the hitting times are approximated
##  with linear interpolation.
			lht.append(
				T[ t - 1 ] + ( T[ t ] - T[ t - 1 ] ) * (
					( j * scale - X[ t - 1 ] ) / ( X[ t ] - X[ t - 1 ] ) ) )
			last_hit = j
## return the crossing information
	if skipFirst :
		lht.pop( 0 )
		lhp.pop( 0 )
	return ( lht, lhp )

## Get the number of subcrossing of the given hiting times and levels.
def get_subxings( ht, hp ) :
	subx = list( )
	prev = ht[ 0 ]
	for l in xrange( 1, len( ht ) ) :
		curr = ht[ l ]
		subxing = list( )
		if curr :
## Find the starting index of the sub-crossings of the current
##  layer l.
			j0 = 0
			while prev[ j0 ] < curr[ 0 ] : j0 += 1
## For each hit time of the crossings of the current layer 
			for i in xrange( 1, len( curr ) ) :
## Find the first index of the next set of sub-crossings.
				j1 = j0 + 1
				while prev[ j1 ] < curr[ i ] : j1 += 1
				subxing.append( j1 - j0 )
## Advance to the set of sub-crossings, which actually belong to
##  the next crossing
				j0 = j1
		prev = curr
		subx.append( subxing )
	return subx

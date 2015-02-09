# -*- coding: UTF-8 -*-
import numpy as np

## The old preamble
## times, points and subcrossings for continuous process
## I. Nazarov & G. Decrouez 2015 -- completely rewrtitten and optimized for extensive simulations
##  X, T -- aligned value and time series of a continuous process
##  X, T = ( {X(t_i)}_{i=1}^n, {t_i}_{i=1}^n )
## X : signal to be analysed
## T : time ( mandatory )
## delta: the base scale of the crossing grid to reconstruct
## deleteFirst: if 1 then delete first crossing at each level
# 
# w: crossing duration ( cell structure, size of the number of levels )
# subx : number of subcrossings ( cell structure, size of the number of levels, first empty )
# hit_point: cell struct., hitting levels 
# hit_time: cell struct., hitting times
# Y.Shen 12.2002, O.D.Jones 7.2003, 6.2007
# P - O. Amblard & G. Decrouez 2011
#
#  [w, subx, hit_point, hit_time] = f_get_w( fname, levels, delta, deleteFirst )


## Compute the crossings of the grid spawned by integers by a process X_T.
def xtree_integer_crossings( T, X ) :
## Assume the process starts at zero, and that its first
##  crossing is at the zero-th line of the grid.
	last_hit = 0
	lht = list([ 0.0 ]) ; lhp = list([ last_hit ])
## Preemptively round the process values up and down
##  in separate arrays to minimize the overhead.
	X_floor = np.floor( X, np.empty_like( X, np.int ) )
	for t in xrange( len( X ) - 1 ) :
		if X[ t ] < X[ t + 1 ] : 	## Upcrossing
## If it was an upcrossing find the smallest integer, larger
##  than the starting point, and the largest integer, smaller
##  than the final point.
			direction = 1
			level0 = X_floor[ t ] + 1 		# np.ceil( X[ t ] )
			level1 = X_floor[ t + 1 ] 		# np.floor( X[ t + 1 ] )
		elif X[ t ] > X[ t + 1 ] :	## Downcrossing
## For a downcrossing find the closest integer less than
##  the origin, and the least integer, greater than the
##  final point.
			direction = -1
			level0 = X_floor[ t ] 			# np.floor( X[ t ] )
			level1 = X_floor[ t + 1 ] + 1 	# np.ceil( X[ t + 1 ] )
		else :
			continue
## It is possible to regard hitting times as a collection of successive
##  stopping times until the first leave of the current grid band:
##    T^n_{k+1} \defn \inf\{ t > T^n_k\,: \, \abs{ X_t - X_{T^n_k} } \geq \delta 2^{-n} \}
##  where T^n_0 = 0.
## It is not a crossing if the process returned to the level
##  of the previous crossing. Adjust the starting level (level0)
		if level0 == last_hit : level0 += direction
## Usually these integers are ordered, but if the process didn't
##  move enough to cross a grid line, that is the endpoints are
##  between two levels but within one band, then these integers
##  would be inverted, when corrected for the direction.
		if ( level0 - level1 ) * direction > 0 : continue
		# beta = ( T[ t + 1 ] - T[ t ] ) / ( X[ t + 1 ] - X[ t ] )
		for level in xrange( level0, level1 + direction, direction ) :
## Due to discretized nature of the process the hitting times are
##  approximated with linear interpolation. The speed is not improved
##  by precalculating the coefficients of the linear interpolation.
			# lht.append( T[ t ] + beta * ( level - X[ t ] ) )
			lht.append( T[ t ] + ( T[ t + 1 ] - T[ t ] ) / ( X[ t + 1 ] - X[ t ] ) * ( level - X[ t ] ) )
## Record the crossed grids: a continuous process during an upcrossing
##  will have almost surely visited all the intermediate grid lines.
			lhp.append( level )
## The last hit is always the last grid line to have been crossed
		last_hit = level1
## Return the crossing times and grid lines
	return ( np.array( lht, np.float ), np.array( lhp, np.int ), np.array( [ ], np.int ) )

## Using the crossing data construct the associated super-crossings
def xtree_super_crossing( T, P, band_width ) :
## By construction, the first crossing is always the zero-th.
	last_hit = 0 ; next_hit = 1
	lht = list( ) ; lhp = list( ) ; subx = list( )
## Instead of initializing the lists with the first hit, add it
	lhp.append( P[ last_hit ] )
	lht.append( T[ last_hit ] )
	while next_hit < len( P ) :
## Find the first time, when the crossing left the ±2 band. This logic
##  critically depends on the assumption that the grid levels are
##  crossed in succession. # Caching the last hit does not help.
		if abs( P[ last_hit ] - P[ next_hit ] ) >= band_width :
## The series of hits between "last_hit" and "next_hit" are
##  the crossings to be aggregated into a super-crossing.
## If we encountered a super-crossing while within the array,
##  commit it to our queue. Note that we do not loose a "crossing"
##  if we just leave this loop on overshooting the length of
##  the array. This is because in this case the last crossing
##  is an incomplete one, and a complete one, with the last within
##  the array, will have been caught.
			lhp.append( P[ next_hit ] )
			lht.append( T[ next_hit ] )
## Count the number of excursions (children in the crossing tree).
			subx.append( next_hit - last_hit )
## Start a new super-crossing
			last_hit = next_hit
		next_hit += 1
	return ( np.array( lht, np.float ), np.array( lhp, np.int ), np.array( subx, np.int ) )

## Adaptive selection of the basic grid scale is based on
##  the standard deviation of the process increments.
## For the given size of the grid band, construct a ten layers deep
##  crossing tree, starting from the root, which corresponds to the
##  coarsest grid.
def xtree_build( T, X, delta = None, max_height = float( 'inf' ) ) :
## By default, delta, the maximum grid spacing, is the standard
##  deviation of the increments.
	delta = np.std( np.diff( X ) ) if delta is None else delta
## Set up the crossing tree structure
	hp = list( ) ; ht = list( ) ; hx = list( )
## Rescale the sample path, so that the grid base scale is 1.0
	Z = ( X - X[ 0 ] ) / delta
## First compute the crossing times and points of the finest
##  integer grid
	lht, lhp, lhx = xtree_integer_crossings( T, Z )
## Add the times and property translated point to the master queue
	hp.append( lhp * delta + X[ 0 ] )
	ht.append( lht )
	hx.append( lhx )
## If the height restriction permits and the crossings did occur
##  iteratively construct crossings of increasingly coarser grids.
	height = 1
	while len( lhp ) > 1 and height < 2**max_height:
		height *= 2
		lht, lhp, lhx = xtree_super_crossing( lht, lhp, height )
		hp.append( lhp * delta + X[ 0 ] )
		ht.append( lht )
		hx.append( lhx )
	return ( ht, hp, hx )

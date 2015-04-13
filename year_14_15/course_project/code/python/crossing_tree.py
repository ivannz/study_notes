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
	X_ceil = np.ceil( X, np.empty_like( X, np.int ) )
	for t in xrange( len( X ) - 1 ) :
		if X[ t ] < X[ t + 1 ] : 	# # Upcrossing
## If it was an upcrossing find the smallest integer, larger
##  than the starting point, and the largest integer, smaller
##  than the final point.
			direction = 1
			level0 = X_ceil[ t ] 			# np.ceil( X[ t ] )
			level1 = X_floor[ t + 1 ] 		# np.floor( X[ t + 1 ] )
		elif X[ t ] > X[ t + 1 ] :	# # Downcrossing
## For a downcrossing find the closest integer less than
##  the origin, and the least integer, greater than the
##  final point.
			direction = -1
			level0 = X_floor[ t ] 			# np.floor( X[ t ] )
			level1 = X_ceil[ t + 1 ] 		# np.ceil( X[ t + 1 ] )
		else :
			continue
## It is possible to regard hitting times as a collection of successive
##  stopping times until the first leave of the current grid band:
##    T^n_{k+1} \defn \inf\{ t > T^n_k\,: \, \abs{ X_t - X_{T^n_k} } \geq \delta 2^{-n} \}
##  where T^n_0 = 0.
## Usually these integers are ordered, but if the process didn't
##  move enough to cross a grid line, that is the endpoints are
##  between two levels but within one band, then these integers
##  would be inverted, when corrected for the direction.
		# if ( level0 - level1 ) * direction > 0 : continue
		# beta = ( T[ t + 1 ] - T[ t ] ) / ( X[ t + 1 ] - X[ t ] )
		for level in xrange( level0, level1 + direction, direction ) :
## It is not a crossing if the process returned to the level
##  of the previous crossing. Adjust the starting level (level0)
			if level != last_hit :
## The last hit is always the last grid line to have been crossed
				last_hit = level
## Due to discretized nature of the process the hitting times are
##  approximated with linear interpolation. The speed is not improved
##  by precalculating the coefficients of the linear interpolation.
			# lht.append( T[ t ] + beta * ( level - X[ t ] ) )
				lht.append( T[ t ] + ( T[ t + 1 ] - T[ t ] ) / ( X[ t + 1 ] - X[ t ] ) * ( level - X[ t ] ) )
## Record the crossed grids: a continuous process during an upcrossing
##  will have almost surely visited all the intermediate grid lines.
				lhp.append( level )
## Return the crossing times and grid lines
	return ( np.array( lht, np.float ), np.array( lhp, np.int ),
		np.array( [ ], np.int ), np.empty( (0,3), np.int ) )

## Using the crossing data construct the associated super-crossings
def xtree_super_crossing( T, P, band_width ) :
## By construction, the first crossing is always the zero-th.
	last_hit = 0 ; next_hit = 1
	lht = list( ) ; lhp = list( ) ; subx = list( ) ; excur = list( )
## Instead of initializing the lists with the first hit, add it
	lhp.append( P[ last_hit ] )
	lht.append( T[ last_hit ] )
	directions = np.sign( np.diff( P ) )
	while next_hit < len( P ) :
## Find the first time, when the crossing left the +\-2 band. This logic
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
## Owen Daffyd Jones 2005: "The first apparent crossing at each level
##  should be excluded, since for a non-Markov process the path from
##  T^n_0 to T^n_1 is not a true crossing" [Citation needed]
			# if last_hit > 0 :
## Count the number of excursions (children in the crossing tree).
			subx.append( next_hit - last_hit )
## Aggregate the directions of excursions: \/ and /\
			# directions = np.sign( np.diff( P[ last_hit : ( next_hit + 1 ) ] ) )
			dirs = directions[ last_hit : next_hit ]
## Count the number of /\ (+-) excluding the last pair //(++) or \\(--)
##  this counting does not depend on the direction of the parent crossing:
##  just +- for U and -+ for the rest, except for -- or ++
			U = np.sum( dirs[ 1:-2:2 ] == -1 )
## The number of \/ is equal to the total number of subcrossings, without
##  the last upward or downward movement.
			excur.append( ( U, ( next_hit - last_hit ) // 2 - U - 1, dirs[ -1 ] ) )
## The format is (/\, \/, ±1) where sign of the last depends on the direction
##  of final crossing.
## Start a new super-crossing
			last_hit = next_hit
		next_hit += 1
	return ( np.array( lht, np.float ), np.array( lhp, np.int ),
		np.array( subx, np.int ), np.array( excur, np.int ) )

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
	hp = list( ) ; ht = list( ) ; hx = list( ) ; ex = list( )
## Rescale the sample path, so that the grid base scale is 1.0
	Z = ( X - X[ 0 ] ) / delta
## First compute the crossing times and points of the finest
##  integer grid
	lht, lhp, lhx, lex = xtree_integer_crossings_superfast( T, Z )
## Add the times and property translated point to the master queue
	hp.append( lhp * delta + X[ 0 ] )
	ht.append( lht )
	hx.append( lhx )
	ex.append( lex )
## If the height restriction permits and the crossings did occur
##  iteratively construct crossings of increasingly coarser grids.
	height = 1
	while len( lhp ) > 1 and height < 2**max_height:
		height *= 2
		lht, lhp, lhx, lex = xtree_super_crossing( lht, lhp, height )
		hp.append( lhp * delta + X[ 0 ] )
		ht.append( lht )
		hx.append( lhx )
		ex.append( lex )
	return ( ht, hp, hx, ex )

def xtree_integer_crossings_fast( T, X ) :
## Assume the process starts at zero, and that its first crossing is at
##  the zero-th line of the grid.
	last_hit = 0
	lht = list([ 0.0 ]) ; lhp = list([ last_hit ])
## Compute the crossing directions
	X_delta = np.diff( X )
	cross_d = np.sign( X_delta, np.empty_like( X_delta, np.int8 ) )
	del X_delta
## Preemptively round the process values down to minimize the overhead.
	X_floor = np.floor( X, np.empty_like( X, np.int ) )
## Note that floor(x) + 1 = ceil(x) for non-integer x only, and for
##  integers by definition floor(x) = ceil(x).
	X_ceil = np.ceil( X, np.empty_like( X, np.int ) )
## If it was an upcrossing find the smallest integer, larger than the starting
##  point, and the largest integer, smaller than the final point. For a
##  downcrossing find the closest integer less than the origin, and the least
##  integer, greater than the final point.
## X_t is rounded either up or down depending on the direction of the crossing
##  it starts: upcrossing -- ceiling, downcrossing -- floor.
	X_begin = np.where( cross_d > 0, X_ceil[:-1], X_floor[:-1] )
## X_{t+1} is rounded down if the crossing it finishes if upward, and is
##  rounded up if the crossing is down, i.e in the direction opposite to the
##  crossing. However due to the nature of xrange() it is also necessary
##  to adjust the final grid level in the direction of the crossing.
	X_final = np.where( cross_d > 0, X_floor[1:], X_ceil[1:] ) + cross_d
	del X_floor, X_ceil
## It is possible to regard hitting times as a collection of successive
##  stopping times until the first leave of the current grid band:
##    T^n_{k+1} \defn \inf\{ t > T^n_k\,: \, \abs{ X_t - X_{T^n_k} } \geq \delta 2^{-n} \}
##  where T^n_0 = 0.
## Within-band excursions and sideways movements are ignored.
	mask = np.logical_and( cross_d != 0, X_begin != X_final )
## Using a X_final[t-1]== X_final[t] rule to eliminate excursions is bad,
##  since sequential small within-band movements may add up to a crossing.
	tau = np.arange( len( X ) - 1, dtype = np.int )[ mask ]
	del mask
	for t in tau :
## Usually X_begin and X_final are ordered or reversed with respect to the
##  crossing direction, but if the process didn't move enough to cross a grid
##  line, i.e. that is the endpoints are between two levels but within one
##  band, then these integers would be inverted, when corrected for the
##  direction.
		for level in xrange( X_begin[ t ], X_final[ t ], cross_d[ t ] ) :
## Either level0 equals the last_hit or no grid line between X_begin and
##  X_final does. 
			if level == last_hit :
## It is not a crossing if the process returned to the level of the previous
##  crossing.
				continue
## Due to discretized nature of the process the hitting times are
##  approximated with linear interpolation. The speed is not improved
##  by precalculating the coefficients of the linear interpolation.
			lht.append( T[ t ] + ( T[ t + 1 ] - T[ t ] ) * ( ( level - X[ t ] ) / ( X[ t + 1 ] - X[ t ] ) ) )
			# lht.append( T[ t ] + ( T[ t + 1 ] - T[ t ] ) * ( ( level - X[ t ] ) / X_delta[ t ] ) )
## Record the crossed grids: a continuous process during an upcrossing
##  will have almost surely visited all the intermediate grid lines.
			lhp.append( level )
## The last hit is always the last grid line to have been crossed
			last_hit = level
## Return the crossing times and grid lines
	return ( np.array( lht, np.float ), np.array( lhp, np.int ),
		np.array( [ ], np.int ), np.empty( (0,3), np.int ) )

def xtree_integer_crossings_superfast( T, X ) :
## Get the closest integers around each value of the process. Note that by definition
##  floor(x) = ceil(x) for integers and floor(x) + 1 = ceil(x) for non-integer x.
	X_floor, X_ceil = np.floor( X ), np.ceil( X )
## Compute movement directions
	cross_d = np.sign( np.diff( X ) )
## For each crossing (X_t, X_{t+1}) determine the intervals potentially crossed
##  by it. For an upcrossing round X_t up (ceiling) and X_{t+1} down (floor),
##  whereas, the levels traversed by a downcrossing are the floor and the ceiling
##  respectively.
	X_begin = np.where( cross_d > 0, X_ceil[:-1], X_floor[:-1] )
	X_end = np.where( cross_d < 0, X_ceil[+1:], X_floor[+1:] )
	del X_floor, X_ceil
## Adjusting the final level in the direction of the crossing helps in detection of
##  the within-band sideways movements, which are non-crossings.
	X_final = X_end + cross_d
## Compute the current sizes of traversals. Later on this data would be used to remove
##  non-crossings and adjust the intially traversed levels for consecutive crossings
##  of the same integer level.
	X_size = np.abs( X_final - X_begin, np.empty_like( cross_d, dtype = np.int ) )
## Eliminate non-movements since their integer bounds and sizes are incorrect.
	X_size[ cross_d == 0 ] = 0
## Within-band sideways movements occur when no integer levels are crossed during
##  the traversal. This happens when floor(x_0) = ceil(x_1)-1 for x_0>x_1, or
##  ceil(x_0) = floor(x_1)+1 for x_0<x_1. In such cases the computed crossing size
##  is zero. The following line masks all movements which do not cross anything.
	is_crossing = X_size > 0
## Find out the indices of movements, which crossed a level. Note, that this array
##  indexes the "is_crossing" array.
	prev = np.concatenate( ( [ 0 ], np.cumsum( is_crossing[ :-1 ], dtype = np.int ) ) ) - 1
## A true crossing of an integer level occurs when this level was not the last level
##  traversed by the crossing preceding this one. This happens if the unadjusted final
##  level of the preceding crossing is equal to the first traversed level of the
##  current one.
	is_recrossing = is_crossing & ( X_end[ is_crossing ][ prev ] == X_begin ) & ( prev >= 0 )
## Effectively, this detects movements, which recrossed the last level of their
##  preceding movement. In order to account for such movements it is necessary to
##  adjust the level each traversal starts with in the direction of the crossing.
## Consecutive re-crossings are properly accounted for, since 
## Adjust initially crossed levels of re-crossings, and update the sizes
	X_begin[ is_recrossing ] += cross_d[ is_recrossing ]
## Adjusting the starting level makes the crossing shorter by one level.
	X_size[ is_recrossing ] -= 1
## Pruning: eliminate one-level traversals which do not contribute a crossing.
	tau = np.arange( len( cross_d ), dtype = np.int )[ X_size > 0 ]
	del prev, is_recrossing
## Arcane programming : get the crossing sizes and compute the index of the end of
##  each crossing in the output array.
	X_size = X_size[ tau ]
	X_index = np.cumsum( X_size, dtype = np.int64 )
## Initialize the final arrays
	X_times, X_values = np.zeros( X_index[ -1 ], np.float ), np.zeros( X_index[ -1 ], np.float )
## Get the index of the beginning of each crossing (group [X_index : X_index + X_size]).
	X_index -= X_size
## In general it is much more likely that a crossing traverses just one or two levels.
##  Invoking the np.arange() for such short events is too expensive. Thus it is more 
##  efficient to handle such cases separately and let numpy's core do all the heavy
## looping. The following selects the times and starting indices of such short crossings.
	tau_first, inx_first = tau[ X_size < 3 ], X_index[ X_size < 3 ]
	tau_second, inx_second = tau[ X_size == 2 ], X_index[ X_size == 2 ]
## Since the process is sampled at discrete moments which are not necessarily spaced
##  uniformly apart, the level hitting times are approximated by linear interpolation.
##  The speed is not improved by the usage of pre-calculated constants. First record
##  the first levels crossed and their times.
	X_values[ inx_first ] = X_begin[ tau_first ]
	X_times[ inx_first ] = T[ tau_first ] + ( T[ tau_first + 1 ] - T[ tau_first ] ) * ( ( X_begin[ tau_first ] - X[ tau_first ] ) / ( X[ tau_first + 1 ] - X[ tau_first ] ) )
	X_values[ inx_second + 1 ] = X_end[ tau_second ]
	X_times[ inx_second + 1 ] = T[ tau_second ] + ( T[ tau_second + 1 ] - T[ tau_second ] ) * ( ( X_end[ tau_second ] - X[ tau_second ] ) / ( X[ tau_second + 1 ] - X[ tau_second ] ) )
	del X_end, tau_first, inx_first, tau_second, inx_second
## Now handle multi-level crossings
	tau, inx, size = tau[ X_size > 2 ], X_index[ X_size > 2 ], X_size[ X_size > 2 ]
	for t, i, j in zip( tau, inx, inx + size ) :
## In case three or more levels were crossed, create an array and use vectorized arithmetic
		lines = np.arange( X_begin[ t ], X_final[ t ], cross_d[ t ], dtype = np.float )
## Add the levels, traversed by the crossing to the final array.
		X_values[ i : j ] = lines
		X_times[ i : j ] = T[ t ] + ( T[ t + 1 ] - T[ t ] ) * ( ( lines - X[ t ] ) / ( X[ t + 1 ] - X[ t ] ) )
	return ( X_times, X_values, np.empty( 0, np.int ), np.empty( (0,3), np.int ) )

## Adaptive selection of the basic (coarsest) grid scale is based on the standard
##  deviation of increments of the sample path of process.
def xtree_build_superfast( T, X, delta = None, max_height = float( 'inf' ) ) :
## Set up the crossing tree structure
	hp = list( ) ; ht = list( ); hx = list( ) ; ex = list( )
## By default, delta, the maximum grid spacing, is the standard
##  deviation of the increments.
	delta = np.std( np.diff( X ) ) if delta is None else delta
## Rescale the sample path, so that the grid base scale is 1.0
	Z = ( X - X[ 0 ] ) / delta
## First compute the crossing times and points of the finest
##  integer grid
	lht0, lhp0, lhx0, lex0 = xtree_integer_crossings_superfast( T, Z )
## Add the times and property translated point to the master queue
	hp.append( lhp0 * delta + X[ 0 ] )
	ht.append( lht0 )
	hx.append( lhx0 )
	ex.append( lex0 )
## If the height restriction permits and the crossings did occur
##  iteratively construct crossings of increasingly coarser grids.
	height = 0
	while len( lhp0 ) > 1 and height < max_height :
## Advnace to the next level of the crossing tree and reduce the scale of
##  the sample path.
		height += 1
		delta *= 2
		Z = ( X - X[ 0 ] ) / delta
## Compute the crossing times and points
		lht1, lhp1, lhx1, lex1 = xtree_integer_crossings_superfast( T, Z )
## Since the crossing times are linearly interpolated, the twofold descaling
##  does not affect the crossing times of the even leves. Thus it is possible
##  to match exactly the times on the consecutive levels of tree.
		if False :
## The very last crossing of the current level has either no nontrvival offspring
##  or is incomplete. 
## If we encountered a super-crossing while within the array,
##  commit it to our queue. Note that we do not loose a "crossing"
##  if we just leave this loop on overshooting the length of
##  the array. This is because in this case the last crossing
##  is an incomplete one, and a complete one, with the last within
##  the array, will have been caught.

			hits = np.searchsorted( lht0, np.concatenate( ( lht1, [ np.infty ] ) ) ) #
			directions = np.concatenate( ( np.sign( np.diff( lhp0 ) ), [ 0.0 ] ) )
		else :
			hits = np.searchsorted( lht0, lht1 )
			directions = np.sign( np.diff( lhp0 ) )
## Owen Daffyd Jones 2005: "The first apparent crossing at each level
##  should be excluded, since for a non-Markov process the path from
##  T^n_0 to T^n_1 is not a true crossing" [Citation needed]
		lht0, lhp0 = lht1, lhp1
		lhit, fhit = hits[+1:]//2 - 1, hits[:-1]//2
		fdir = directions[ 1::2 ]
## Aggregate the directions of up-down excursions /\. The down-up |/ are computed
##  based on these and the number of pairs of sub-crossings.
		mask = ( fdir == -1 )
		mask[ lhit ] = False
		ud_excur = np.cumsum( mask )
## The format is (/\, \/, ±1) where sign of the last depends on the direction
##  of final crossing.
		excursions = np.empty( ( len( hits ) - 1, 3 ), np.int )
## Count the number of /\ (+-) excluding the last pair //(++) or \\(--)
##  this counting does not depend on the direction of the parent crossing:
##  just +- for U and -+ for the rest, except for -- or ++
		excursions[:,0] = ud_excur[ lhit ] - np.where( fhit > 0, ud_excur[ fhit-1 ], 0 )
## The number of \/ is equal to the total number of subcrossings, without
##  the last upward or downward movement.
		excursions[:,1] = lhit - fhit - excursions[:,0]
## Record the direction
		excursions[:,2] = fdir[ lhit ]
		hp.append( lhp0 * delta + X[ 0 ] )
		ht.append( lht0 )
## Count the number of offspring of the current scale crossings. This is just
##  the number of subcrossings between two consecutive lower-scale crossings. 
		hx.append( np.array( np.diff( hits ), np.int ) )
		ex.append( excursions )
	return ( ht, hp, hx, ex )

################################################################################
################################# CODE REVIEW! #################################
################################################################################
### Checking the original code by Decrouez and Owen:
### Results identical!! Checked through monte-carlo assertions of identity.
def f_get_w_int( T, X, deleteFirst = False ) :
#    if deleteFirst
	if deleteFirst :
#        last_hit = 0;
		last_hit = 0
#	else
	else :
#        last_hit = 1;
		last_hit = 1
#    end
# h_t=zeros(1,lx*ceil(max(diff(y)))); %% at this scale, upper bound on the total # of crossings
# h_p=h_t;
	lht = list( ) ; lhp = list( )
# compt=1;
# y_floor = floor( y );
	X_floor = np.floor( X, np.empty_like( X, np.int ) )
	X_ceil = np.ceil( X, np.empty_like( X, np.int ) )
# for i = 1:(lx-1)
	for t in xrange( len( X ) - 1 ) :
#     if y(i) ~= y(i+1)
		if X[ t ] == X[ t + 1 ] :
			continue
#         if y(i) < y(i+1)
		if X[ t ] < X[ t + 1 ] :
#             step = 1;
			direction = 1
#             x_init = ceil(y(i));
			level0 = X_ceil[ t ]
#             x_final = floor(y(i+1));
			level1 = X_floor[ t + 1 ]
#         else
		elif X[ t ] > X[ t + 1 ] :
#             step = -1;
			direction = -1
#             x_init = floor(y(i));
			level0 = X_floor[ t ]
#             x_final = ceil(y(i+1));
			level1 = X_ceil[ t + 1 ]
#         end
#         for j = x_init:step:x_final
		# has_run = False
		for level in xrange( level0, level1 + direction, direction ) :
			# has_run = True
#             if j ~= last_hit
			if level != last_hit :
#                 h_t(compt) = t(i) + (j - y(i))/(y(i+1) - y(i))*(t(i+1) - t(i));
				lht.append( T[ t ] + ( level - X[ t ] ) / ( X[ t + 1 ] - X[ t ] ) * ( T[ t + 1 ] - T[ t ] ) )
#                 h_p(compt) = j*delta*scale + x(1);
				lhp.append( level )
#                 compt=compt+1;             
				last_hit = level
#                 last_hit = j;
#             end
		# if has_run :
		# 	assert( level == level1 )
#         end
#     end
# end
	return ( np.array( lht, np.float ), np.array( lhp, np.int ),
		np.array( [ ], np.int ), np.empty( (0,3), np.int ) )

def f_get_w( T, X, levels = [ ], delta = 1.0, deleteFirst = False ) :
# hit_point=cell(length(levels),1);
# hit_time=cell(length(levels),1);
# w=hit_time;
# comptscale=1;
## Set up the crossing tree structure
	hp = list( ) ; ht = list( )
# z = ( x - x( 1 ) ) / delta ;
	Y = ( X - X[ 0 ] ) / delta
	for n in levels :
# for n = levels
#     scale = 2^n;
		scale = float( 2**n )
#     y = z/scale;
		Z = Y / scale
# .....
		lht, lhp, lhx = f_get_w_int( T, Z, deleteFirst )
#     hit_point{comptscale}=h_p(1:compt-1);
		hp.append( lhp * scale * delta + X[ 0 ] )
#     hit_time{comptscale}=h_t(1:compt-1);
		ht.append( lht )
#     w{comptscale}=diff(h_t(1:compt-1));
#     comptscale=comptscale+1;   
# hit0=[hit_time{1}' hit_point{1}'];
	hit0 = ht[ 0 ]
# subx=cell(length(levels),1);
	hx = list( [np.array( [], np.int )] )
# for level = 2:length(levels)   
	for level in xrange( 1, len( levels ) ) :
#    hit1=[hit_time{level}' hit_point{level}'];  
		hit1 = ht[ level ]
#    if ~isempty(hit1)
		if len( hit1 ) > 0 :
#        j0 = 1;
			j0 = 0
#        sx=zeros(1,size(hit1,1)-1);
			sx = list()
#        compt=1;
#        while hit0(j0,1) ~= hit1(1,1), j0 = j0 + 1; end
			while ( j0 < len( hit0 ) ) and ( hit0[ j0 ] != hit1[ 0 ] ) : j0 += 1
#        for i = 2:size(hit1,1)
			for i in xrange( 1, len( hit1 ) ) :
#            j1 = j0 + 1;
				j1 = j0 + 1
#            while hit0(j1,1) ~= hit1(i,1), j1 = j1 + 1; end           
				while ( j1 < len( hit0 ) ) and ( hit0[ j1 ] != hit1[ i ] ) : j1 += 1
				#if j1 < len( hit0 ) :
#            sx(compt)=j1-j0;
				sx.append( j1 - j0 )
#            compt=compt+1;
#            j0 = j1;
				j0 = j1
#        end
#	end
#    subx{level}=sx;
		hx.append( np.array( sx, np.int ) )
#    hit0 = hit1;
		hit0 = hit1
	return ( ht, hp, hx )
#end

		






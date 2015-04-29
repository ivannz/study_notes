# -*- coding: UTF-8 -*-
import numpy as np

def xtree_integer_crossings( T, X ) :
## Compute movement directions
	X_direction = np.sign( np.diff( X ) )
## Get the closest integers around each value of the process. Note that by definition
##  floor(x) = ceil(x) for integers and floor(x) + 1 = ceil(x) for non-integer x.
	X_floor, X_ceil = np.floor( X ), np.ceil( X )
## For each crossing (X_t, X_{t+1}) determine the intervals potentially crossed
##  by it. For an upcrossing round X_t up (ceiling) and X_{t+1} down (floor),
##  whereas, the levels traversed by a downcrossing are the floor and the ceiling
##  respectively.
	X_begin = np.where( X_direction > 0, X_ceil[:-1], X_floor[:-1] )
	X_final = np.where( X_direction < 0, X_ceil[+1:], X_floor[+1:] )
	del X_floor, X_ceil
## Compute the current sizes of traversals. Later on this data would be used to remove
##  non-crossings and adjust the intially traversed levels for consecutive crossings
##  of the same integer level.
## Adjusting the final level in the direction of the crossing helps in detection of
##  the within-band sideways movements, which are non-crossings.
	size = np.abs( X_final + X_direction - X_begin, np.empty_like( X_direction, dtype = np.int ) )
## Eliminate non-movements since their integer bounds and sizes are incorrect.
	size[ X_direction == 0 ] = 0
## Within-band sideways movements occur when no integer levels are crossed during
##  the traversal. This happens when floor(x_0) = ceil(x_1)-1 for x_0>x_1, or
##  ceil(x_0) = floor(x_1)+1 for x_0<x_1. In such cases the computed crossing size
##  is zero. The following line masks all movements which do not cross anything.
	is_crossing = size > 0
## Find out the indices of movements, which crossed a level. Note, that this array
##  indexes the "is_crossing" array.
	previous_crossing = np.concatenate( ( [ 0 ], np.cumsum( is_crossing[ :-1 ], dtype = np.int ) ) ) - 1
## A true crossing of an integer level occurs when this level was not the last level
##  traversed by the crossing preceding this one. This happens if the unadjusted final
##  level of the preceding crossing is equal to the first traversed level of the
##  current one.
	is_recrossing = is_crossing & ( X_begin == X_final[ is_crossing ][ previous_crossing ] ) & ( previous_crossing >= 0 )
## Effectively, this detects movements, which recrossed the last level of their
##  preceding movement. In order to account for such movements it is necessary to
##  adjust the level each traversal starts with in the direction of the crossing.
####  Consecutive re-crossings are properly accounted for, since ...
## Adjust initially crossed levels of re-crossings, and update the sizes.
	X_begin[ is_recrossing ] += X_direction[ is_recrossing ]
## Adjusting the starting level necessarily makes the crossing shorter by one level.
	size[ is_recrossing ] -= 1
	del previous_crossing, is_recrossing
## Pruning: eliminate one-level traversals which do not contribute a crossing because
##  they are recrossings.
	tau = np.nonzero( size > 0 )[ 0 ] # np.arange( len( X_direction ), dtype = np.int )[ size > 0 ]
## Get the sizes of genuine crossings
	size = size[ tau ]
## Arcane programming : compute the index of the end of each crossing in the output
##  array. Hopefully the incre,emts are not so large as to cause 64-bit int overflow.
	inx = np.cumsum( size, dtype = np.int64 )
## Initialize final crossing times and values arrays
	X_times, X_values = np.zeros( inx[ -1 ], np.float ), np.zeros( inx[ -1 ], np.float )
## Get the index of the beginning of each crossing (aka group [inx:inx+size] in
##  python colon notation).
	inx -= size
## In general it is much more likely that a crossing traverses just one or two levels.
##  Invoking the np.arange() for such short events is too expensive. Thus it is
##  more efficient to handle such cases separately and let numpy's core do all the
##  heavy looping. The following selects the times and starting indices of such
##  short crossings. 
	tau1, inx1 = tau[ size < 3 ], inx[ size < 3 ]
	tau2, inx2 = tau[ size == 2 ], inx[ size == 2 ] + 1
## Since the process is sampled at discrete moments which are not necessarily spaced
##  uniformly apart, the level hitting times are approximated by linear interpolation.
##  The speed is not improved by the usage of pre-calculated constants. First record
##  the first levels crossed and their times.
	X_values[ inx1 ] = X_begin[ tau1 ]
	X_times[ inx1 ] = T[ tau1 ] + ( T[ tau1 + 1 ] - T[ tau1 ] ) * ( ( X_begin[ tau1 ] - X[ tau1 ] ) / ( X[ tau1 + 1 ] - X[ tau1 ] ) )
	X_values[ inx2 ] = X_final[ tau2 ]
	X_times[ inx2 ] = T[ tau2 ] + ( T[ tau2 + 1 ] - T[ tau2 ] ) * ( ( X_final[ tau2 ] - X[ tau2 ] ) / ( X[ tau2 + 1 ] - X[ tau2 ] ) )
	del tau1, inx1, tau2, inx2
## Now handle multi-level crossings
	tau, inx, size = tau[ size > 2 ], inx[ size > 2 ], size[ size > 2 ]
	for t, i, j in zip( tau, inx, inx + size ) :
## In case three or more levels were crossed, create an array and use vectorized arithmetic
		lines = np.arange( X_begin[ t ], X_final[ t ] + X_direction[ t ], X_direction[ t ], dtype = np.float )
## Add the levels, traversed by the crossing to the final array.
		X_values[ i : j ] = lines
		X_times[ i : j ] = T[ t ] + ( T[ t + 1 ] - T[ t ] ) * ( ( lines - X[ t ] ) / ( X[ t + 1 ] - X[ t ] ) )
	return X_times, X_values

## Adaptive selection of the basic (coarsest) grid scale is based on the standard
##  deviation of increments of the sample path of process.
def xtree_build( T, X, delta = None, max_height = float( 'inf' ) ) :
## Set up the crossing tree structure
	hp = list( ) ; ht = list( ) ; hx = list( ) ; ex = list( ) ; wt = list( )
## By default, delta, the maximum grid spacing, is the standard
##  deviation of the increments.
	delta = np.std( np.diff( X ) ) if delta is None else delta
## Rescale the sample path, so that the grid base scale is 1.0
	Z = ( X - X[ 0 ] ) / delta
## First compute the crossing times and points of the finest
##  integer grid
	lht0, lhp0 = xtree_integer_crossings( T, Z )
## Add the times and property translated point to the master queue
	ht.append( lht0 )
	hp.append( lhp0 * delta + X[ 0 ] )
	hx.append( np.empty( 0, np.int ) )
	ex.append( np.empty( ( 0, 3 ), np.int ) )
	wt.append( np.empty( 0, np.int ) )
## If the height restriction permits and the crossings did occur
##  iteratively construct crossings of increasingly coarser grids.
	height = 0
	while len( lhp0 ) > 1 and height < max_height :
## Advance to the next level of the crossing tree and reduce the scale of the sample
##  path.
		height += 1 ; delta *= 2
		Z = ( X - X[ 0 ] ) / delta
## Compute the crossing times and points
		lht1, lhp1 = xtree_integer_crossings( T, Z )
## Owen Daffyd Jones 2005: "The first apparent crossing at each level should be
##  excluded, since for a non-Markov process the path from T^n_0 to T^n_1 is not
##  a true crossing" [Citation needed]
		# ToDo: Implement this option!
## Since the crossing times are linearly interpolated, the twofold descaling does not
##  affect the crossing times of the even levels. Thus it is possible to match exactly
##  the times on the consecutive levels of tree. The following logic depends on the
##  condition that lht0[0] <= lht1[0]. In fact theoretically lht1 is a subset of lht0.
		hits = np.searchsorted( lht0, lht1 ) # np.concatenate( ( lht1, [ np.infty ] ) )
## Find out the index of the last pair of subcrossings before the current crossing
		last_hit = hits // 2 - 1
## Compute the direction of the second subcrossing in each pair
		directions = np.sign( np.diff( lhp0 ) )[ 1::2 ]
## Locate true up-down excursions (down-down)
		up_down_mask = directions < 0
		up_down_mask[ last_hit[ 1: ] ] = False
## Count the number of up-down excursions encountered so far.
		up_down_total = np.concatenate( ( [ 0 ], np.cumsum( up_down_mask )[ last_hit[ 1: ] ] ) )
		up_down_total[ last_hit < 0 ] = 0
## The very last crossing of the current level is either degenerate (no offspring)
##  or is incomplete. An incomplete crossing is one with known beginning but unknown
##  end. Such last crossing occurs only an the end of the examined sample path and
##  is likely to include a sub-crossing which might not even be paired!
## The format is [#/\, #\/, ±1] where sign of the last depends on the direction of
##  the final crossing. The number of subcrossing on the current level is one less
##  than the number of hitting times of the current level (coincides with len(hits)).
## Aggregate the directions of up-down excursions /\. The down-up |/ are computed
##  based on these and the number of pairs of sub-crossings given by lhit - fhit.
		excursions = np.empty( ( len( hits ) - 1, 3 ), np.int )
## Count the number of up-down /\ excluding the last pair up-up // or down-down \\.
## The number of down-up \/ is equal to the total number of excursions, without
##  the up-down ones /\.
		excursions[:,0] = np.diff( up_down_total )
		excursions[:,1] = np.diff( hits ) // 2 - excursions[:,0] - 1
		excursions[:,2] = directions[ last_hit[ 1: ] ]
## Commit the current level to the queue
		ht.append( lht1 )
		hp.append( lhp1 * delta + X[ 0 ] )
## Count the number of offspring of the current scale crossings. This is just
##  the number of subcrossings between two consecutive lower-scale crossings. 
		hx.append( np.array( np.diff( hits ), np.int ) )
		ex.append( excursions )
		wt.append( np.diff( lht1 ) )
		lht0, lhp0 = lht1, lhp1
	return ( ht, hp, hx, ex, wt )

# a = list( )
# ud = 0 ; du = 0
# directions = np.sign( np.diff( lhp0 ) )
# for f, s in zip( directions[0::2], directions[1::2] ) :
#     if f > 0 and s < 0 :
#         ud += 1
#     elif f < 0 and s > 0 :
#         du += 1
#     else :
#         a.append( ( ud, du, s ) )
#         ud = du = 0
# aa = np.array( a )
# print np.allclose( aa, excursions )


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
# ## [j,j+i,j+2i, ...,j+m*i] j:i:k and m = fix((k-j)/i)
# ## i=1 : [j,j+1,j+2, ...,k-1,k]
# ## i=-1: [j,j-1,j-2, ...,k+1,k]
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
	hp = list( ) ; ht = list( ) ; hv = list( )
# z = ( x - x( 1 ) ) / delta ;
	Y = ( X - X[ 0 ] ) / delta
	for n in levels :
# for n = levels
#     scale = 2^n;
		scale = float( 2**n )
#     y = z/scale;
		Z = Y / scale
# .....
		lht, lhp, lhx, lhv = f_get_w_int( T, Z, deleteFirst )
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
	return ( ht, hp, hx, hv )
#end


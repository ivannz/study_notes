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
## Owen Daffyd Jones 2005: "The first apparent crosing at each level
##  should be excluded, since for a non-Markov process the path from
##  T^n_0 to T^n_1 is not a true crossing" [Citation needed]
			# if last_hit > 0 :
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
	lht, lhp, lhx = xtree_integer_crossings_fast( T, Z )
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

def xtree_integer_crossings_fast( T, X ) :
## Assume the process starts at zero, and that its first crossing is at
##  the zero-th line of the grid.
	last_hit = 0
	lht = list([ 0.0 ]) ; lhp = list([ last_hit ])
## Compute the crossing directions
	X_delta = np.diff( X )
	cross_d = np.sign( X_delta, np.empty_like( X_delta, np.int ) )
## Preemptively round the process values down to minimize the overhead.
	X_floor = np.floor( X, np.empty_like( X, np.int ) )
## Note that floor(x) + 1 = ceil(x) for non-integer x only, and for
##  integers by definitions floor(x) = ceil(x).
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
	for t in xrange( len( X ) - 1 ) :
		if cross_d[ t ] == 0 : continue
## Usually X_begin and X_final are ordered or reversed with respec to the
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
			lht.append( T[ t ] + ( T[ t + 1 ] - T[ t ] ) * ( ( level - X[ t ] ) / X_delta[ t ] ) )
## Record the crossed grids: a continuous process during an upcrossing
##  will have almost surely visited all the intermediate grid lines.
			lhp.append( level )
## The last hit is always the last grid line to have been crossed
			last_hit = level
## Return the crossing times and grid lines
	return ( np.array( lht, np.float ), np.array( lhp, np.int ), np.array( [ ], np.int ) )

################################################################################
################################# CODE REVIEW! #################################
################################################################################
### Checking the original code by Decrouez and Owen:
### Results identical!! Checked through monet-carlo assertions of indenity.
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
	return ( np.array( lht, np.float ), np.array( lhp, np.int ), np.array( [ ], np.int ) )

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

		






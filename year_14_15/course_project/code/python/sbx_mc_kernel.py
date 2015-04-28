import numpy as np
from scipy.stats import chisqprob as cdchisq

from fgn import fbm
from xtree import xtree_build

generator = fbm( N = 2**21+1, H = .65 )
generator.set_rnd( np.random.RandomState( 123 ) )
T, X = generator( )

max_levels, max_crossings = 10, 10

delta = np.std( np.diff( X ) )

Tnk, hp, Znk, ex, Wt = xtree_build( T, X, delta = delta )
# Tnk, hp, Znk, ex, Wt = xtree_build( T, X, delta = delta, max_height = max_levels + 1 )

Nn = np.array( [len( Tk ) for Tk in Tnk ], dtype = np.int ).reshape( ( len( Tnk ), -1 ) )
counts = np.zeros( ( max_levels + 1, max_crossings // 2 ), dtype = np.float )
for level, subcrossings in enumerate( Znk[ 1: ], 0 ) :
	sx_count, sx_freq = np.unique( subcrossings, return_counts = True )
	sx_count = np.minimum( sx_count, max_crossings )
	level = np.minimum( level, max_levels )
	counts[ level, sx_count // 2 - 1 ] = sx_freq

levels = range( 1, np.minimum( 4, len( Nn ) ) )
p_count = np.sum( counts[ levels,: ], axis = 0 )
p_distr = p_count / p_count.sum( )
distr = counts / counts.sum( axis = 1 ).reshape( counts.shape[ 0 ], -1 )
stat = np.sum( Nn[ levels ] * ( ( distr[ levels,: ] - p_distr )**2 / p_distr ) )
dof = len( levels ) * ( distr.shape[ 1 ] - 1 )
cdchisq( stat, dof )



####################################################################################################
## The Monte Carlo kernel performs a single experiment and returns its
## K'xT slice of the result.
def monte_carlo_kernel( generator, K, T, delta = None ) :
	from scipy.stats import chisqprob as cdchisq
## K is the maximal number of subcrossings beyond which the data
##  on subcrossings is agregateed into the tail. (truncation parameter)
	K = 2 * ( K // 2 + 1 )
## T is the number of detailed levels of the crossing tree to stored in
##  output. Any crossings of grids coarser than the treshold are
##  aggregated.
	distr = np.zeros( ( T+1, K // 2 ), dtype = np.int )
## Generate a replication of the process
	t, x = generator( )
## G. Decrouez 2015-02-12: the selection of the spacing of the finest grid
##  based on the scale of the process is crucial as it allows comparison
##  of the crossing tree between different sample paths.
	if delta is None :
# 		sigma = np.std( np.diff( x ) )
## Set the base scale so that the constructed tree is 6 layers deep.
# 		delta = 2**( np.floor( np.log2( np.max( x ) - np.min( x ) ) ) - 6 )
		delta = None
## 2015-04-08 : I think 6-level deep crossing tree is too poor for any analysis
##  which is why it is necessary to take the variance into account
## Get the hitting times and points
	ht, hp, subx, excur = xtree_build_superfast( t, x, delta )
## Pool all the crossing counts together and construct the distribution matrix 
## ToDo
## ChiSquare test for distributional match: choose the levels to test the hypothesis on
	levels = range( -5, -3+1, 1 )
## The total number of crossings of \\delta 2^n grid
	Nn = np.array( [ len( ht[ l ] ) for l in levels ] ).reshape( (len(levels), -1 ) )
## The number of subcrossings of size \\delta 2^{n-1} that
##  make up the k-th crossing of n-th level.
	Zn = [ subx[ l ] for l in levels ]
## Compute the empirical distribution of the number of subcrossings in the pooled sample
	bins, f_pool = np.unique( np.concatenate( tuple( Zn ) ), return_counts = True )
	p_pool = f_pool / np.sum( f_pool, dtype = np.float64 )
## Calculate the empirical probabilites at each level
	p_lvl = np.zeros( ( len( Nn ), len( bins ) ), dtype = np.float )
	for l, z in enumerate( Zn, 0 ) :
		cat, freq = np.unique( z, return_counts = True )
		p_lvl[ l, np.searchsorted( bins, cat ) ] = freq / np.sum( freq, dtype = np.float64 )
## Run the chi_square test on each level of the tree
	stat = np.sum( Nn * ( ( p_lvl - p_pool )**2 / p_pool ) )
	dof = len( levels ) * ( len( bins ) - 1 )
## Now summarize the distribution: Simplify
## ToDo
	for level, xing in enumerate( subx[1:], 0 ) :
## Count the absolute frequencies
		c, f = np.unique( xing, return_counts = True )
## We collect the data on the distribution truncated by K (including).
## Thus truncate the number of subcrossings by K -- everything less is
##  detailed, otherwise -- agregated into tail
		c = np.minimum( c, K )
## Similarly truncate the level by T + 1
		if level >= T : level = T
## Fill the offspring distribution matrix
		for c, f in zip( c, f ) :
			distr[ level, c // 2 - 1 ] += f
## Use straightforward bivariate regression implementation.
## We could also use ML estimators, which is much better!
## See the handwritten notes.
## Estimate the conditional distribution of excursions
## lvl x (++ --) x (+- -+) # (/\\, \\/, +/-1) ==  (u,d,s)
	hatvee = np.zeros( (len( excur ), 2, 2 ), dtype = np.int )
	for i, lex in enumerate( excur, 0 ) :
		if lex.shape[0] < 1 : continue
		hatvee[i,0,:] = lex[ lex[:,2]>0, :2 ].sum( axis = 0 )
		hatvee[i,1,:] = lex[ lex[:,2]<0, :2 ].sum( axis = 0 )
## Estimate the probability distribution of the number of offspring
## ToDo
	return distr, (stat,dof, cdchisq(stat, dof)), hatvee


####################################################################################################
op = { 'L' : 10, 'K' : 10 }

## max_levels -- the number of level of the tree to construct (from the finest grid)
##  beyond this level (L+1, L+2, ...) the levels are pooled.
## max_crossing -- the threshold for the number of subcrossings beyond which they
##  are considered to be in the tail. The distribution data would contain bins:
##  {2}, {4}, ... {2k}, ... {2K, 2(K+1), ...}, where K = max_crossing
max_levels, max_crossings = op.get( 'L', 6 ), op.get( 'K', 20 )

## Tnk[n][k] -- the time (approximate) of the k-th crossing of the gird with resolution
##  \delta 2^n. XTnk[n][k] -- the value of the process at the time of the crossing:
##  alwayts equal to the shifted and scaled grid level crossed.

## Znk[n][k] -- the number of subcrossings of a finer grid (crossings of
##  \delta 2^{n-1}) that make up the k-th corssing of a coarser grid (\delta 2^n).
##  Undefined (empty array) for n = 0.

## Vnk[n][k] -- the number of up-down and down-up excursions in the k-th
##  crossing of size \delta 2^n. Undefined for n = 0.

## Wnk[n][k] -- the waiting time between the k-th and k+1-st crossing of the grid
##  of spacing \delta 2^n
Tnk, XTnk, Znk, Vnk, Wnk = xtree_build( T, X, delta = delta )

## Get the total number of crosssings of \delta 2^n resolution
## Nn[n] -- the total number of crossings of grid with spacing \delta 2^n
Nn = np.array( [ len( Tk ) - 1 for Tk in Tnk ],
	dtype = np.float ).reshape( ( len( Tnk ), -1 ) )

## Dnm[n][m] -- m<M : the total number of crossings of grid \delta 2^{n+1} with 2(m+1)
##  subcrossings of grid \delta 2^n. The values in column M are the number
##  of crossings of \delta 2^{n+1} with not less than 2(M+1) subcrossings.
## 32bit Integers should be enough
Dnm = np.zeros( ( max_levels + 1, max_crossings // 2 ), dtype = np.float )
for n, Zk in enumerate( Znk[ 1: ], 0 ) :
	n = max_levels if n > max_levels else n
	Z_count, Z_freq = np.unique( Zk, return_counts = True )
## Truncate the observed crossings
	Z_count = np.minimum( Z_count, max_crossings )
	mask = ( Z_count < max_crossings )
	Dnm[ n, Z_count[ mask ] // 2 - 1 ] += Z_freq[ mask ]
	Dnm[ n, max_crossings // 2 - 1 ] += np.sum( Z_freq[ ~mask ] )

## VDn[n][d][e] -- the total number of up-down(e=0) and down-up(e=1) excursions
##  (/\ and \/ subcrosssings od \delta 2^n respectively) in a downward(d=0) or
##  upward(d=1) crossing of spacing \delta 2^{n+1} (level n+1). Levels beyond
##  max_levels are agregated.
VDn = np.zeros( ( max_levels + 1, 2, 2 ), dtype = np.int )
for n, Vk in enumerate( Vnk[ 1: ], 0 ) :
	n = max_levels if n > max_levels else n
	VDn[ n, 0 ] += np.sum( Vk[ Vk[ :, 2 ] < 0 ], axis = 0 )[:2]
	VDn[ n, 1 ] += np.sum( Vk[ Vk[ :, 2 ] > 0 ], axis = 0 )[:2]

## Consdier the first H bins
num_bins = 2
levels = range( 4, 8 )

Z_hist = np.zeros( ( len( levels ), num_bins+1 ), dtype = np.float )
Z_hist[:, :num_bins ] = Dnm[ levels, :num_bins ]
Z_hist[:, num_bins ] = np.sum( Dnm[ levels, num_bins: ], axis = 1 )

P_pooled = np.sum( Z_hist, axis = 0 ) / np.sum( Z_hist )
P_emp = Z_hist / Nn[1:][levels]

chi_stat = np.sum( Nn[1:][ levels ] * ( P_emp - P_pooled )**2 / P_pooled )
chi_dof = ( P_emp.shape[0] - 1 ) * ( P_emp.shape[1] - 1 )

####################################################################################################
def mc_kernel( generate_sample, **op ) :
## Get the desired resolution
	max_levels, max_crossings = op.get( 'L', 6 ), op.get( 'K', 20 )
## Generate a replication of the process -- a sample path
	T, X = generate_sample( )
## G. Decrouez 2015-02-12: the selection of the spacing of the finest grid based on
##  the scale of the process is crucial as it allows comparison of the crossing tree
##  between different sample paths.
	delta = np.std( np.diff( X ) )
## 2015-04-08 : I think 6-level deep crossing tree is too poor for any analysis which
##  is why it is necessary to take the variance, i.e the inherent scale, into account.
## Generate a replication of the process -- a sample path
	T, X = generate_sample( )
## G. Decrouez 2015-02-12: the selection of the spacing of the finest grid based on
##  the scale of the process is crucial as it allows comparison of the crossing tree
##  between different sample paths.
	delta = np.std( np.diff( X ) )
## 2015-04-08 : I think 6-level deep crossing tree is too poor for any analysis which
##  is why it is necessary to take the variance, i.e the inherent scale, into account.
	Tnk, hp, Znk, ex, Wt = xtree_build( T, X, delta = delta, max_height = max_levels + 1 )
## The total number of crossings of \delta 2^n grid
	Nn = np.array( [len( Tk ) for Tk in Tnk ], dtype = np.int ).reshape( ( len( Tnk ), -1 ) )
## Summarize the distribution: 32-bit integers should be enough 
	counts = np.zeros( ( max_levels + 1, max_crossings // 2 ), dtype = np.float )
	for level, subcrossings in enumerate( Znk[ 1: ], 0 ) :
## Compute the frequencies of subcrossings
		sx_count, sx_freq = np.unique( subcrossings, return_counts = True )
## Put the largest numbers into one bin
		sx_count = np.minimum( sx_count, max_crossings )
## Agregate the top (least populated) levels of the tree.
		level = np.minimum( level, max_levels )
		counts[ level, sx_count // 2 - 1 ] = sx_freq
## Compute the empirical distribution of the number of subcrossings in the pooled sample
	levels = range( 1, np.minimum( 4, len( Nn ) ) )
	p_count = np.sum( counts[ levels,: ], axis = 0 )
	p_distr = p_count / p_count.sum( )
	distr = counts / counts.sum( axis = 1 ).reshape( counts.shape[ 0 ], -1 )
	stat = np.sum( Nn[ levels ] * ( ( distr[ levels,: ] - p_distr )**2 / p_distr ) )
	dof = len( levels ) * ( distr.shape[ 1 ] - 1 )
	return counts, ( stat, dof )

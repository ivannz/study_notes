# -*- coding: UTF-8 -*-
import numpy as np
from scipy.stats import chisqprob as cdchisq

## The procedure below runs a Chi-square distribution test on the specified levels
##  of the summarized crossing tree to see whether the relevant subcrossings have
##  the same distribution. Parameters "p" and "q" specify the first and the last
##  level respectively, and "h" sets the number of bins to runs the Chi-square test on.
def pooled( Djnk, Njn, H = None, p = 5, q = 9, h = 5 ) :
## Set the array of layers
	levels = np.arange( p, 1 + q, dtype = np.int ) - 1
## Extract the required levels of the tree: compute the per level normalisation
##  constant and the empirical distribution.
	Cjn = np.reshape( np.sum( Djnk[ :, levels ], axis = ( 2, ), dtype = np.float ), ( Djnk.shape[ 0 ], q - p + 1, -1 ) )
	phat_jnk = Djnk[ :, levels ] / Cjn
## Get the necessary layers and compute the pooled subcrossing distribution.
##  In the cube Djnk the first layer of $\delta 2^{-1}$ subcrossings is omitted:
##  Djnk[:,0,:] are the frequencies of subcrossings of size $\delta 2^0$ in a single
##  crossing of lower resolution resolution $\delta 2^{0+1}$.
	pool_jk = np.sum( Djnk[ :, levels ], axis = ( 1, ) ) / np.sum( Cjn, axis = ( 1, ) )
## Re-bin the data: per level empirical distribution	
	phat = np.zeros( phat_jnk.shape[ :-1 ] + ( h + 1, ), dtype = np.float )
	phat[:,:,:-1], phat[:,:,-1] = phat_jnk[:,:,:h], np.sum( phat_jnk[:,:,h:], axis = ( 2, ) )
## Re-bin the data: empirical distribution of pooled levels	
	pool = np.zeros( ( Djnk.shape[ 0 ], 1, h + 1 ), dtype = np.float )
	pool[:,0,:-1], pool[:,0,-1] = pool_jk[:,:h], np.sum( pool_jk[:,h:], axis = ( 1, ) )
## Compute the chi-squared test-statistic for each replication. However the array
##  of the number of crossings of grid \delta 2^n has the n=0 crossings.
	stat = np.sum( Njn[ :, levels + 1 ] * ( ( phat - pool )**2 ) / pool, axis = ( 1, 2, ) )
	dof  = np.array( phat.shape[ 0 ] * [ ( q - p ) * ( h - 1 ) ], dtype = np.int )
## Return the result of the Chi-square pooling test for each Monte-Carlo replication.
	return stat, dof, cdchisq( stat, dof )

## The procedure below test whether the subcrossing distribution of the pooled tree
##  levels (from p to q inclusive) significantly diverges from the conjectured one.
def theoretical( Djnk, Njn, H, p = 5, q = 9, h = 5 ) :
	levels = np.arange( p, 1 + q ) - 1
## Pool the specified levels of the crossing tree to test the conjectured hypothesis
	Cj = np.reshape( np.sum( Djnk[ :, levels ], axis = ( 1, 2, ), dtype = np.float ), ( Djnk.shape[ 0 ], -1 ) )
	Nj = np.sum( Njn[ :, levels + 1 ], axis = 1 )
	pool_jk = np.sum( Djnk[ :, levels ], axis = ( 1, ) ) / Cj
## Re-bin the distribution
	pool = np.zeros( ( Djnk.shape[ 0 ], h + 1 ), dtype = np.float )
	pool[:,:-1], pool[:,-1] = pool_jk[:,:h], np.sum( pool_jk[:,h:], axis = ( 1, ) )
## Compute the theoretical probabilities
	theta = 2.0 ** ( 1.0 - 1.0 / H )
	prob = np.concatenate( ( theta * ( 1 - theta ) ** np.arange( h ), [ ( 1 - theta ) ** h ] ) )
## Compute the \chi^2 test statistic
	stat = np.sum( Nj * ( ( pool - prob )**2 ) / prob, axis = ( 1, ) )
	dof  = np.array( pool.shape[ 0 ] * [ h - 1 ], dtype = np.int )
## Return the result of the Chi-square theoretical distribution test for each
##   Monte-Carlo replication.
	return stat, dof, cdchisq( stat, dof )

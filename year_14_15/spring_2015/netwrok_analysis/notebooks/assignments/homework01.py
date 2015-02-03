import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as la

## Construct a regression model
def lm_model( X, Y, intercept = True ) :
    T = np.array( Y, dtype = float )
    M = np.array( X, dtype = float )
    if intercept is True :
        M = np.vstack( [ np.ones( len( Y ) ), M ] ).T
	return (M,T, intercept)

## Define the OLS regression routine:
def lm_fit( model ) :
	M, T, intercept = model
	MMinv = la.inv( ## implement (X'X)^{-1} (X'Y)
		np.dot( M.T, M ) ) 
	coef = np.dot( MMinv,
		np.dot( M.T, T ) )
## Estimate the residual standard deviation
	resid = T - np.dot(M, coef)
	dof = len( T ) - len( coef )
	RSS = np.dot( resid.T, resid )
	return (coef, RSS, dof, MMinv )

#####################################################################
#+ 0. Load the data (yes, it is a milestone!)
## Load the word count dataset
wordcount = np.fromregex(
	'./data/wordcounts.txt', r"(\d+)\s+(.{,32})",
	[ ( 'freq', np.int64 ), ( 'word', 'S32' ) ] )


#####################################################################
##+ 1. Check that Zipf's Law holds
## Pre-sort the frequencies: in ascending order of frequencies
wordcount.sort( order = 'freq' )
freqs = wordcount[ 'freq' ]
## PRoduce ranks: from 1 up to |W|
ranks = np.arange( 1, len( wordcount ) + 1, dtype = float )[::-1]
## The probability of a word frequency being not less than the 
##  frequency of a gien word w it exactly the ratio of the w's rank
##  to the total number of words.
probs = ranks / len( wordcount )

## estimate f_k\sim C k^{-\gamma} model
mdl = lm_model( np.log( ranks ), np.log( freqs ), True )
coef, rss, dof, XX = lm_fit( mdl )

## Define the fitted Zipf's law
# zipf = lambda r : np.exp( coef.dot( ( 1, np.log( r ) ) ) )
zipf = lambda r : np.exp( coef[0] + coef[1] * np.log( r ) )

## Show how well is was estimated.
plt.loglog( freqs, probs, "xr" )
plt.plot( zipf( ranks ), probs, "-b" )
plt.xlabel( "frequency" ) ; plt.ylabel( "ranks" )
plt.title( "Wordcount data" )
plt.show( )

######################################################################
##+ 2. Assuming that the data is distributed according to the Power Law, find
##  * $\alpha$ of the distribution
##  * mean sample variance $\sigma^2$

## ML estimator of the power law in the "tail" (x≥u):
##  x_k \sim C x^{-\alpha} 1_{[u,+∞)}(x).
def mle_alpha( data, threshold ) :
## Keep the data observations, that we consider to be in the tail
	tail = np.array( [ v for v in data if v >= threshold ] )
## Estimate the mean log of the peaks over threshold
	sum_log = np.sum( np.log( tail ) ) / ( len( tail ) + 0.0 )
## Use the closed form expression for the value of the power at an optimum
	alpha = 1.0 + 1.0 / ( sum_log - np.log( threshold ) )
## Using the delta-method compute the s.e of the estimate.
	return alpha, ( alpha - 1 ) / np.sqrt( len( tail ) )

## Get the ML estimate
alpha_ml, alpha_ml_sd = mle_alpha( freqs, freqs.min( ) )

## Let's suppose that the rank is proportional to the complementary CDF
##  of a power law: $\bar{F}(x) = {\left(\frac{x}{u}\right)}^{1-\alpha}$
##  Thus the following econometric model is to be estimated:
##  $\log \text{rank} \sim C + (1-\alpha) \log \text{freq} + \epsilon$
mdl = lm_model( np.log( freqs ), np.log( ranks ), True )
beta, rss, dof, XX = lm_fit( mdl )
## Transform the coefficient
alpha_ls = 1 - beta[ 1 ]

## The regression estimate of the power should be close
##  to the ML estimate
print "the OLS estimate of alpha is %f\n" % alpha_ls
print "Whereas the ML estimate is %f (%f) \n" % ( alpha_ml, alpha_ml_sd )
print "Since ML is more theoretically sound, the relative error is %f%%\n" % (
	100 * np.abs( 1.0 - alpha_ls / alpha_ml ), )

## The mean and the sample variance of the sample
##  frequency distribution:
print "The average frequency over the sample is ", freqs.mean(), "\n"
print "The sample variance is ", freqs.var(), "\n"

## Theoretical mean and variance of the power law distribution
##  significantly depend on the power parameter.
## Indeed for $x\sim \frac{\alpha-1}{u} {\left( \frac{x}{u} \right)}^\alpha$ one has the following:
##   $E(x) = \frac{\alpha-1}{\alpha-2} u$ if $\alpha>2$
##   $E(x^2) = \frac{\alpha-1}{\alpha-3} u^2$ if $\alpha>3$
## The estimated parameter is less than 2, implying that the frequency
##  distribution is unlikely to have even a finite mean under the
##  assumed distribution.

#####################################################################
##+ 3. Produce summary of the frequencies: min, max, mean, median
## Does it make sense to compute these summaries? What does the mean frequency tell us?
print "The minimum frequency is ", freqs.min(), "\n"
print "The mean frequency is ", freqs.mean(), "\n"
print "The median frequency is ", np.median( freqs ), "\n"
print "The maximum frequency is ", freqs.max(), "\n"

######################################################################
## A necessary preamble: function definitions.

## Define a convenience function for estimating the power parameter
##  of the continuous power law
from scipy.stats import kstest
def ks_dist( data, threshold ) :
## Estimate the power given the current threshold
	alpha, sd = mle_alpha( data, threshold )
## Construct the CDF in the current environment
	cdf = lambda x : 1.0 - ( x / threshold ) ** ( 1.0 - alpha )
## Return the output of the out-of-the box Kolmogorov-Smirnov test:
##  the infinity norm of the difference between the distribution functions.
	d, pv = kstest( [ v for v in data if v >= threshold ], cdf )
	return (d, pv), (alpha, sd)

## Two functions below implement the same functionality as the previous two
##  but instead of the continuous version they work with the discrete power law.
from scipy.special import zeta
from scipy.optimize import minimize
## The discrete power law gives marginally different results
##  \Pr(N=n) \defn \frac{1}{\zeta(\gamma)} n^{-\gamma}, n -- positive integer
def mle_alpha_d( data, threshold ) :
## Keep the data observations, that we consider to be in the tail
	tail = np.array( [ v for v in data if v >= threshold ] )
## Estimate the mean log of the peaks over threshold
	sum_log = np.sum( np.log( tail ) ) / ( len( tail ) + 0.0 )
## Define minus log-likelihood of the discrete power law
	loglik = lambda alpha : np.log( zeta( alpha, threshold ) ) + alpha * sum_log
## Compute the ML estimate of the exponent, with a view to using it as the
##  initial seed for the numerical minimizer for better convergence.
	alpha_0 = 1.0 + 1.0 / ( sum_log - np.log( threshold ) )
	res = minimize( loglik, ( alpha_0, ), method = 'Nelder-Mead', options = { 'disp': False } )
## Return the "optimal" argument, regardless of its quality. Potentially DANGEROUS!
	return res.x[ 0 ], float( 'nan' )

def ks_dist_d( data, threshold ) :
## Estimate the power given the current threshold
	alpha, sd = mle_alpha_d( data, threshold )
## Construct the CDF in the current environment
	cdf = lambda k : zeta( alpha, threshold + k ) / zeta( alpha, threshold )
## Return the output of the out-of-the box Kolmogorov-Smirnov test:
##  the infinity norm of the difference between the distribution functions.
	d, pv = kstest( [ v for v in data if v >= threshold ], cdf )
	return (d, pv), (alpha, sd)


# In my opinion, mentioning the KS test itself, not just it's statistic, is fraught with danger of misuse. Indeed, in some young minds it might plant an idea that it is the possibile to use the P-value produced to test the hypothesis of goodness-of-fit. The distribution of the KS statistic under the null hypothesis rests heavily on the assumption of absolute continuity of the hypothesised distribution. 
# by test when comparing between the models with different thresholds. This is 
# The  KS test and  seems a little
# The main issue with 

#####################################################################
## These helper functions invert an array and count the number of
##  occurrences of distinct values in an array
def values( data, frequency = False ) :
	bins = dict( )
## For each value in the given array, add the index of each occurrence
##  into the bin dedicated to the encountered value.
	for i, x in enumerate( sorted( data ) ) :
## Prepend the current occurrence of a value, unless it has never been
##  seen before, in which case initialise the list of indices for it.
		bins[ x ] = bins.get( x, [] ) + [ i ]
	return bins

def counts( data ) :
## Count the number of times a value occurs in the array.
	counts = dict( )
	for x in data :
## If the values has not been seen yet, then initialize it to
##  a single occurrence otherwise increment its counter.
		counts[ x ] = counts.get( x, 0 ) + 1
	return counts.items( )

def mean_excess( data ) :
	data = np.array( sorted( data, reverse = True ) )
## Compute the last positions in the sorted array of each repeated observation
	ranks = rankdata( data, method = 'max' )
## Since the array is sorted, the number of observation exceeding the current
##  is givne by difference between the length of the array and the max-rank.
	excesses = np.array( np.unique( len( data ) - ranks ), dtype = np.int )
## Get the thresholds
	thresholds = data[ excesses ]
## Get the sum of all values greater than the current threshold 
	mean_excess = np.cumsum( data )[ excesses ] / ( excesses + 0.0 ) - thresholds
	return np.array( thresholds, mean_excess )

#####################################################################
## + 0. Read the graph
## Load the network routing graph first as it is the smallest. It is
##  an undirected graph.
import networkx as nx
# G = nx.read_edgelist( "./data/fb_Princeton.txt", create_using = nx.Graph( ) );
G = nx.read_edgelist( "./data/web_Stanford.txt", create_using = nx.Graph( ) );
# G = nx.read_edgelist( "./data/network.txt", create_using = nx.Graph( ) );
node_degree = G.degree( )
deg = np.array( node_degree.values( ), dtype = np.int )

#####################################################################
##+ 1. Are they correspondent to power law?

## First let's draw the frequency plot of the node degree distribution.
degree_freq = np.array( counts( deg ) )

plt.title( "Node degree frequency" )
plt.xlabel( "degree" ) ; plt.ylabel( "frequency" )
plt.loglog( degree_freq[:,0], degree_freq[:,1], "rx" )
plt.show( )

## The empirical degree distribution does not correspond to a power
##  law per se, but it definitely has some heavy tailed behaviour,
##  which exhibits itself, when the only data exceeding same truncated
##  is considered.

## To see this, let's plot the complimentary cumulative distribution
##  function for various tail thresholds.
def ccdf( data, threshold ) :
## Count the occurrences of values over some threshold in the array
	freq = np.array( counts(
			[ v for v in data if v >= threshold ] ),
		dtype = float )
## Sort the counts along the growing values they correspond to
	freq = freq[ freq[ :, 0 ].argsort( ), : ]
## ... and compute the fraction of data with values lower than the current
	freq[:,1] = 1.0 - np.cumsum( freq[ :,1 ], dtype = float ) / sum( freq[ :,1 ] )
	return freq

cc_10 = ccdf( deg, 10 )
cc_100 = ccdf( deg, 100 )

plt.title( "Truncated degree cCDF" )
plt.xlabel( "degree" ) ; plt.ylabel( "probability" )
plt.loglog( cc_10[:,0], cc_10[:,1], "r-", linewidth = 2 )
plt.loglog( cc_100[:,0], cc_100[:,1], "k-", linewidth = 2 )
plt.show( )

## Clearly the chances of an extremely high node degree decay proportional
##  to the value of the degree on a log-log scale.

#####################################################################
##+ 2. Find max and mean values of incoming and outcoming node degrees
## Since the network graph is undirected it does not make sense to
##  distinguish in- and out- nodes. Thus let's check the range of the
##  general (two-way) degree.

print "The degrees range from %d to %d" % ( min( deg ), max( deg ) ) #, "\n"
print "The average degree over the sample is %.3f" % ( G.size( ) / G.order( ) ) #, "\n"
print "The degree standard deviation is %.3f" % ( np.sqrt( np.var( deg ) ) ) #, "\n"
print "The median degree is %d" % ( np.median( deg ) ) #, "\n"

#####################################################################
##+ 3. Find $\alpha$ via Maximum Likelihood and calculate $\sigma^2$
##+ 4. Determine $x_{min}$ via Kolmogorov-Smirnov test

## We have reasons to believe there are some power law-like effects in
##  the behaviour of the node degree (treated as a random variable).
##  Let's pursue this lead and estimate the exponent in the power law
##  and select the most likely breakpoint, beyond which the degree
##  is heavy tailed.

#####################################################################
## Get the ML estimate of the exponent parameter.
alpha_ml, alpha_ml_se = mle_alpha( deg, min( deg ) )
print "The Maximum likelihood estimate of the exponent of the node degree distribution is %.3f (%.4f)\n" % ( alpha_ml, alpha_ml_se )

#####################################################################
## Run the KS threshold selection routine
thresholds = np.unique( deg )
## The ks_dist() function returns a tuple of the following parameters:
##  * ( KS-distance, PV of the KS-test ), ( MLE of alpha, the standard error of the MLE )
ks_min = np.array( [ ks_dist( deg, u ) for u in thresholds ] )

## Produce the hill plot: the correspondence between the threshold
##  and the estimated exponent.
plt.figure( 1, figsize = ( 5, 5 ) )
plt.label( 'The hill plot of the degree distribution' )
plt.plot( thresholds, ks_min[:,1,0], "r-")

## In fact the KS-metric is the $L^\infty$ norm on the set of distribution
##  functions.
plt.figure( 2, figsize = ( 5, 5 ) )
plt.label( 'The KS metric of the estimated and the empirical distributions' )
plt.plot( thresholds, ks_min[:,0,0], "r-")
plt.show( )

## Select the x_min that brings the KS metric to its minimum on the given
##  degree data. Note the first threshold is removed, since it is likely
##  to yield very biased estimate.
i_min = np.argmin( ks_min[10:,0,0] )
x_min = thresholds[ i_min ]
alpha_ml, alpha_ml_se = ks_min[ i_min, 1, : ]
print "The Kolmogorov-Smirnov metric yielded %.1f as the optimal threshold\n\n" % ( x_min)
print "'Optimal' exponent is %.3f (%.3f)\n\n" % ( alpha_ml, alpha_ml_se )

x = np.arange( x_min, 2 * np.max( deg ) )
deg_ccdf = ccdf( deg, x_min )
pwr_ccdf = lambda x : ( x / ( x_min + 0.0 ) ) ** ( 1.0 - alpha_ml )

plt.title( "Truncated degree cCDF" )
plt.xlabel( "degree" ) ; plt.ylabel( "probability" )
plt.plot( x, pwr_ccdf( x ), "b-", linewidth = 2 )
plt.plot( deg_ccdf[:,0], deg_ccdf[:,1], "r-", linewidth = 2 )
plt.show( )

## network routing -- exhibits Pareto tails
## fb_Princeton, web_Stanford -- exhibit anything but the power law


import numpy as np
import np.random as rnd
import matplotlib.pyplot as plt

def counts( data ) :
## Count the number of times a value occurs in the array.
	counts = dict( )
	for x in data :
## If the values has not been seen yet, then initialize it to
##  a single occurrence otherwise increment its counter.
		counts[ x ] = counts.get( x, 0 ) + 1
	return np.array( counts.items( ) )

theta, alpha = 0.01, 0.25
## Generate geometrically distributed random variate, and then
##  based on each one gemerate a variate from the binomial distribution.
Geom = rnd.geometric( theta, 1000000 )
Thin = [ rnd.binomial( n - 1, alpha ) for n in Geom ]

xx, tt = counts( Geom ), counts( Thin )
plt.plot( tt[:,0], tt[:,1], "b-", linewidth = 2 )
plt.plot( xx[:,0], xx[:,1], "r-", linewidth = 2 )
plt.show()

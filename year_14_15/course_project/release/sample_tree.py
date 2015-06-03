import os, re
import numpy as np
from numpy.random import RandomState

from processes.fgn_numpy import fgn
from processes.fbm import fbm

from crossing_tree import xtree_build
import matplotlib.pyplot as plt

## Set the folder where to dump the plots
plot_path = os.path.realpath( "./output/" )

## Sertup the path ;ength and the Hurst exponent
N, H = 2**18+1, 0.55

gen = fbm( fgn, N, H, time = True )
gen.initialize( RandomState( 123 ) )

## Generate a sample path
T, X = gen( )
# Rescale the process for illustrative purposes

## Set the base scale to the median
delta = np.median( np.abs( np.diff( X ) ) )

## Build a crossing tree
Tnk, Xnk, Znk, Vnk, Wnk = xtree_build( T, X, delta = delta )

l = len( Tnk ) - 2

## Plot the sample path                                      
figure = plt.figure( figsize = ( 12, 9 ) )
axis = figure.add_subplot( 111 )
axis.set_xlim( -0.01, 1.01 )
axis.set_xticks( Tnk[l-3], minor = True )
axis.set_xticks( np.linspace( 0, 1, num = 11 ), minor = False )
axis.set_xticklabels( [ "%0.2f"%(t,) for t in np.linspace( 0, 1, num = 11 ) ], minor = False )
axis.set_yticks( X[ 0 ] + np.arange( -20, 20 ) * delta * 2 ** ( l - 3 ) )
axis.plot( T, X, linestyle = '-', color = 'gray', label = 'X(t)' )
# axis.set_color_cycle( plt.cm.rainbow( np.linspace( 0, 1, 4 )[ ::-1 ] ) )
axis.plot( Tnk[l-3], Xnk[l-3], '-o' )
axis.plot( Tnk[l-2], Xnk[l-2], '-o' )
axis.plot( Tnk[l-1], Xnk[l-1], '-o' )
axis.plot( Tnk[l-0], Xnk[l-0], '-o' )
axis.grid( )
axis.grid( b = True, which = 'minor', color = 'black', linestyle = ':')
# plt.show()
plt.savefig( os.path.join( plot_path, "sample_path.pdf" ) , format = 'pdf' )


## Create the crossing tree
figure = plt.figure( figsize = ( 12, 9 ) )
axis = figure.add_subplot( 111 )
axis.set_xlim( -0.01, 1.01 )
levels = np.arange( l-3, l+1, dtype = np.int )
axis.set_xticks( Tnk[ levels[ 0 ] ], minor = True )
axis.set_xticks( np.linspace( 0, 1, num = 11 ), minor = False )
axis.set_xticklabels( [ "%0.2f"%(t,) for t in np.linspace( 0, 1, num = 11 ) ], minor = False )
axis.set_ylim( 0.9 * delta * 2**( l - 3 ),  1.1 * delta * 2** ( l + 1 ) ) ; axis.set_yscale( 'log' )
axis.set_yticks( delta * 2**levels )
for j in levels :
    lht0, lht1 = Tnk[j], Tnk[j+1]
    parent = np.searchsorted( lht1, lht0, 'left' ) - 1 ; parent[ parent < 0 ] = 0
## Draw the line segments between two levels
    axis.plot( [ lht1[ parent ], lht0 ],
        [ len( lht0 ) * [ delta * 2.0**(j+1) ], len( lht0 )*[ delta * 2.0**j ] ],
        '-o', color = 'black' )

axis.grid( )
# axis.grid( b = True, which = 'minor', color = 'black', linestyle = ':')
axis.set_title( r"""An example of the crossing tree""" )
# plt.show( )
plt.savefig( os.path.join( plot_path, "sample_tree.pdf" ) , format = 'pdf' )



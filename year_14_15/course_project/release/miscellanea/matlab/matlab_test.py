import numpy as np
from numpy.random import RandomState

from processes.fgn_pyfftw import fgn
from processes.fbm import fbm

import scipy.io as io

from main import path_analyse

N, H = 2**20 + 1, 0.55
generator = fbm( fgn, N = N, H = H, time = True )
generator.initialize( RandomState( 123 ) )

T, X = generator( )

io.savemat( './output/data.mat', { 'T': T, 'X': X, 'H': H, 'N': N } )

# import numpy as np
# import scipy.io as io
# mat = io.loadmat( './output/data.mat' )
# H, X, T, N = mat['H'][0], mat['X'][0], mat['T'][0], mat['N'][0]

from crossing_tree import xtree_build, f_get_w
## Matlab's and Numpy's std differ : one of them is biased, whereas the other is not.
delta = 1e-3 # np.std( np.diff( X ) )
Tnk, Xnk, Znk, Vnk, Wnk = xtree_build( T, X, delta = delta )
print np.unique( Znk[1], return_counts = True )

## %% clear all ; format long ;
## %% 
## %% cd c:\study_notes\year_14_15\course_project\code
## %% 
## %% load('C:\study_notes\year_14_15\course_project\code\project\output\data.mat')
## %% 
## %% Zr = 2 : 2 : 40 ;
## %% delta = 1e-3 ; % std( diff( X ) ) ;
## %% [ w, subx, hp, ht ] = f_get_w( X, T, [ 0 : 16 ], delta, 0 ) ;
## %% Z = [ subx{ 1+1 } ] ;
## %% for z = Zr
## %% 	sum( Z == z ) %/ length( Z )
## %% end

if False :
##  # ht, hp, hx, hv = f_get_w( T, X, range( 0, 17 ), delta )
##  # print np.all( [ np.allclose( a0,a1 ) for a0, a1 in zip( Xnk, hp ) ] )
  io.savemat( './output/xtree.mat', { 'Xnk': Xnk, 'Tnk': Tnk } )
##
  Z = ( X - X[ 0 ] ) / delta
  Z_floor = np.floor( Z, np.empty_like( Z, np.float64 ) )
  Z_ceil = np.ceil( Z, np.empty_like( Z, np.float64 ) )
  io.savemat( './output/ceil_floor.mat', {
    'py_ceilz': Z_ceil,
    'py_floorz': Z_floor,
    } )

################################################################################
# delta = np.std( np.diff( X ) )
# delta = np.subtract( *np.percentile( np.diff( X ), [ 75, 25 ] ) )
# delta = np.subtract( *np.percentile( X, [ 100, 0 ] ) ) * 2**( -6 )
delta = np.median( np.abs( np.diff( X ) ) )
Nn, Dnk, Vnde, ( Wnp, Wbar_n, Wstd_n ) = path_analyse(
              T, X, delta = delta, max_levels = 20, max_crossings = 40 )

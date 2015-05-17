import numpy as np
from numpy.random import RandomState

from fgn import fbm
import scipy.io as io

N, H = 2**20 + 1, 0.5
generator = fbm( N = N, H = H, time = True )
generator.set_rnd( RandomState( 123 ) )

T, X = generator( )
io.savemat( './output/data.mat', { 'T': T, 'X': X, 'H': H, 'N': N } )

# import numpy as np
# import scipy.io as io
# mat = io.loadmat( './output/data.mat' )
# H, X, T, N = mat['H'][0], mat['X'][0], mat['T'][0], mat['N'][0]

from crossing_tree import xtree_build, f_get_w
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
delta = np.std( np.diff( X ) )
Tnk, Xnk, Znk, Vnk, Wnk = xtree_build( T, X, delta = delta )

Nn = np.zeros( ( 1 + max_levels + 1, 1 ), dtype = np.int )
for n, Xk in enumerate( Xnk, 0 ) :
  n = max_levels + 1 if n > max_levels + 1 else n
  Nn[ n ] += len( Xk ) - 1

Dnk = np.zeros( ( max_levels + 1, max_crossings // 2 ), dtype = np.int )
for n, Zk in enumerate( Znk[ 1: ], 0 ) :
  n = max_levels if n > max_levels else n
  Z_count, Z_freq = np.unique( Zk, return_counts = True )
  Z_count = np.minimum( Z_count, max_crossings )
  mask = ( Z_count < max_crossings )
  Dnk[ n, Z_count[ mask ] // 2 - 1 ] += Z_freq[ mask ]
  Dnk[ n, max_crossings // 2 - 1 ] += np.sum( Z_freq[ ~mask ] )

Vnde = np.zeros( ( max_levels + 1, 2, 2 ), dtype = np.int )
for n, Vk in enumerate( Vnk[ 1: ], 0 ) :
  n = max_levels if n > max_levels else n
  Vnde[ n, 0 ] += np.sum( Vk[ Vk[ :, 2 ] < 0 ], axis = 0 )[:2]
  Vnde[ n, 1 ] += np.sum( Vk[ Vk[ :, 2 ] > 0 ], axis = 0 )[:2]

prc = np.array( [ 0.5, 1.0, 2.5, 5.0, 10, 25, 50, 75, 90, 95, 97.5, 99, 99.5 ] )
Wnp = np.zeros( ( max_levels, ) + prc.shape, dtype = np.float )
Wbarn = np.zeros( ( max_levels, 1 ), dtype = np.float )
Wstdn = np.zeros( ( max_levels, 1 ), dtype = np.float )
for n, Wk in enumerate( Wnk[1:], 0 ) :
  if len( Wk ) and n < max_levels :
    Wbarn[ n ], Wstdn[ n ], Wnp[ n ] = np.average( Wk ), np.std( Wk ), np.percentile( Wk, prc )

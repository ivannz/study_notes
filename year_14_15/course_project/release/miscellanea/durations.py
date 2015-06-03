
L, H, Njn, Djnk, Vjnde, delta_j, Wjnp, Wbarjn, Wstdjn = results['med']['FBM'][2]

L, H, Njn, Djnk, Vjnde, delta_j, Wjnp, Wbarjn, Wstdjn = results['med']['FBM'][4]



Wnp = np.average( Wjnp, axis = 0 ) * 2**(
	-np.arange( Wjnp.shape[1], dtype = np.float  ).reshape( ( Wjnp.shape[ 1 ], 1 ) ) / H )

## Show a plot for average crossing durations 
figure = plt.figure( figsize = ( 16, 9 ) )
axis = figure.add_subplot( 111 )
axis.set_xlabel( 'Level of the crossing tree' )
axis.set_yscale( 'log' )
axis.set_ylabel( r"""$2^{-n H^{-1}}\mathbb{E}{W}_n$""", rotation = 90 )
axis.boxplot( [ Wbarjn[:,n] * 2 ** ( - n / H ) for n in range( Wbarjn.shape[ 1 ] ) ] )
plt.show( )

# w_n^H x(1) ~ x(w_n) ~ 2^n x(w_1) ~ 2^n w_1^H x(1)
# w_n = (4^{\frac{1}{2H}})^n w_1
p, q = 7, 8
levels = np.arange( p, q + 1 ) - 1

Wnp = np.average( Wjnp, axis = 0 ) * 2**(
	-np.arange( Wjnp.shape[1], dtype = np.float  ).reshape( ( Wjnp.shape[ 1 ], 1 ) ) / H )



scaled_wjnp = Wjnp * 2**(
	-np.arange( Wjnp.shape[1], dtype = np.float  ).reshape( ( 1, Wjnp.shape[ 1 ], 1 ) ) / H )

figure = plt.figure( figsize = ( 16, 9 ) )
axis = figure.add_subplot( 111 )
# axis.set_xticklabels( [ "%0.1f" %(q,) for q in [ 0.5, 1.0, 2.5, 5.0, 10, 25, 50, 75, 90, 95, 97.5, 99, 99.5 ] ], minor = False )
axis.set_xlabel( 'Empirical of crossing durations' )
axis.set_yscale( 'log' )
axis.set_ylabel( r"""Quantile""", rotation = 90 )
axis.boxplot( [ scaled_wjnp[:,p:q+1,k] for k in range( scaled_wjnp.shape[ 2 ] ) ] )
plt.show( )


plt.show()

layers = list( )
for n in range( Wjnp.shape[ 1 ] ) :
	layers.append( Wbarjn[:,n] * 2 ** ( - n / H ) )

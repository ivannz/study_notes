
L, H, Njn, Djnk, Vjnde, Wjnp, Wbarjn, Wstdjn = results['med']['FBM']

L, H, Njn, Djnk, Vjnde, Wjnp, Wbarjn, Wstdjn = results['med']['HRM-3'][2]



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
p, q = 7,7
levels = np.arange( p, q + 1 ) - 1

def figure_09( figure, method, kind, p, q ) :
	figure = plt.figure( figsize = ( 16, 9 ) )
	axis = figure.add_subplot( 111 )
	axis.set_ylabel( 'Probability' )
	axis.set_xlabel( r"""Quantile""", rotation = 00 ) ; axis.set_xscale( 'log' )
	axis.set_color_cycle( plt.cm.rainbow( np.linspace( 0, 1, 5 )[::-1] ) )
	for L, H, Njn, Djnk, Vjnde, Wjnp, Wbarjn, Wstdjn in results[ method ][ kind ] :
		# Wnp = np.average( Wjnp, axis = (0,) ) * 2**(
		# 	-np.arange( Wjnp.shape[1], dtype = np.float  ).reshape( ( Wjnp.shape[ 1 ], 1 ) ) / H )
		# for n in levels :
		# 	axis.plot( Wnp[n], [ 0.5, 1.0, 2.5, 5.0, 10, 25, 50, 75, 90, 95, 97.5, 99, 99.5 ],
		# 		linestyle = '-', marker = 'o', label = L )
		scaled_Wjnp = Wjnp * 2**(
			-np.arange( Wjnp.shape[1], dtype = np.float  ).reshape( ( 1, Wjnp.shape[ 1 ], 1 ) ) / H )
		Wp = np.average( scaled_Wjnp[:,levels], axis = (0, 1, ) )
		axis.plot( Wp, [ 0.5, 1.0, 2.5, 5.0, 10, 25, 50, 75, 90, 95, 97.5, 99, 99.5 ],
			linestyle = '-', marker = 'o', label = L )
	axis.set_title( r"""Empirical distribution of crossing durations at levels %d-%d""" %(p,q,) )
	legend = axis.legend( loc = 'upper left', frameon = 1 )
	frame = legend.get_frame( )
	frame.set_facecolor( 'whitesmoke' )

plt.show( )




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

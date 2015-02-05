# -*- coding: UTF-8 -*-
import numpy as np

## fonction de Weirstrass.
##
## n : nombre de points
## nu0 : frequence rŽduite de base
## nue : freq d'Žchantillonnage
## H : paramtre de rugositŽ
## choix : cha”ne de caractre 
##         'd' : weierstrass dŽterministe
##         'a' : weierstrass alŽatoire
##
## [x,t] = weir( N,H,nu0,nue,choix );
from numpy.random import uniform as unif
def synthweier( N, H, nu0, nue, deterministic = False, seed = None ) :
## Create thee time axis
	t = np.arange( N, dtype = np.float ) / nue
	nmax = int( np.fix( np.log( 0.5 * nue ) / np.log( nu0 ) ) )
## Initialize the trajectory
	x = np.zeros( N, dtype = np.float )
## Generate pahse shifts
	phi = np.zeros( 2 * nmax + 1 )
	if not deterministic :
		if seed is not None : rnd.seed( seed )
		phi = 2 * np.pi * unif( size = 2 * nmax + 1 )
## Construc the Weierstrass process
	for i in xrange( - nmax, nmax + 1 ) :
		x += nu0 ** ( - i * H ) * ( np.cos( phi[ i + nmax ] ) -
			np.cos( phi[ i + nmax ] + 2 * np.pi * ( nu0 ** i ) * t ) )
	return ( x, t )

## Example
# import matplotlib.pyplot as plt
# x, t = synthweier( 2**15, .7, 1.2, 1000 )
# plt.plot(t, x, "b-", linewidth = 2 )
# plt.show( )


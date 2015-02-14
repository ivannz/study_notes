# -*- coding: UTF-8 -*-
import numpy as np

import numpy.random as rnd

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
def synthweier( N, H, nu0, nue, deterministic = False, seed = None ) :
## I don't understand the mathematics behind these prameters -- must
##  investigate further!
## Create the time axis.
	t = np.arange( N, dtype = np.float ) / nue
	nmax = int( np.fix( np.log( 0.5 * nue ) / np.log( nu0 ) ) )
## Initialize the trajectory
	x = np.zeros( N, dtype = np.float )
## Generate pahse shifts
	phi = np.zeros( 2 * nmax + 1 )
	if not deterministic :
		if seed is not None : rnd.seed( seed )
		phi = 2 * np.pi * rnd.uniform( size = 2 * nmax + 1 )
## Construct the Weierstrass process
	for i in xrange( - nmax, nmax + 1 ) :
		x += nu0 ** ( - i * H ) * ( np.cos( phi[ i + nmax ] ) -
			np.cos( phi[ i + nmax ] + 2 * np.pi * ( nu0 ** i ) * t ) )
	return ( t, x )

class synth_Weier(object):
	def __init__(self, N, H):
		super(synth_Weier, self).__init__()
		self.__N = N
		self.__H = H
	def __call__( self, seed = None ) :
		return synthweier( self.__N, self.__H, nu0 = 1.2, nue = self.__N - 1, seed = seed )
	def reset( self ):
		pass

## Example
# import matplotlib.pyplot as plt
# x, t = synthweier( 2**15, .7, 1.2, 1000 )
# plt.plot(t, x, "b-", linewidth = 2 )
# plt.show( )

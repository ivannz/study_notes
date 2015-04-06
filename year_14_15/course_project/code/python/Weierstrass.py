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

class synth_Weier_2(object):
	def __init__(self, N, nu0, H ):
		super(synth_Weier_2, self).__init__()
## Initialize the number and the array of layers
		M = int( np.fix( np.log( N / 2.0 ) / np.log( nu0 ) ) )
		self.__layers = np.arange( - M, M + 1, dtype = np.int )
## Setup cached values
		self.__multo = nu0 ** ( - H * self.__layers )
		self.__multi = 2 * np.pi * ( nu0 ** self.__layers )
## Initalize the time grid [0,1]
		self.__tau = np.linspace( 0, 1, num = N + 1, dtype = np.float)
	def __call__( self, seed = None, deterministic = False ) :
## The implementation below basically constructs a fine mesh on the time axis $(t_k)_{i=0}^N\\in [0,1]$ with  
##   $$0 = t_0< t_1 < \\ldots < t_{N-1} < t_N = 1$$
##  and then evaluates the $W_\\alpha$ at the knots:
##   $$W_i = W_\\alpha(t_i) = \\sum_{k=-M}^{M} \\nu_0^{-k\\alpha} \\bigl( \\cos(\\phi_k) - \\cos(\\nu_0^k 2 \\pi t_i + \\phi_k) \\bigr)$$
##  where $M = \\bigl[\\frac{\\log \\frac{1}{2} \\nu_e}{\\log \\nu_0}\\bigr]$ and
##   $(\\phi_k)_{k=-M}^M\\sim \\mathcal{U}[0,1]$ -- random phase shift of the lay
## This array holds the values of Weierstrass function at every knot of the mesh on [0,1]
		W = np.zeros( len( self.__tau ), dtype = np.float )
## Initialize the phase shifts
		phi = np.zeros( len( self.__layers ), dtype = np.float )
		if not deterministic :
			if seed is not None : rnd.seed( seed )
			phi = 2 * np.pi * rnd.uniform( size = len( self.__layers ) )
## pre-compute the cosines of the phase shift
		cosphi = np.cos( phi )
## For each layer evaluate the trigonometric funcitons
		for m in range( len( self.__layers ) ) :
			w_layer = np.cos( self.__multi[ m ] * self.__tau + phi[ m ] )
			W += self.__multo[ m ] * ( cosphi[ m ] - w_layer )
		return self.__tau, W
	def reset( self ):
		pass

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

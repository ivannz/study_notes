# -*- coding: UTF-8 -*-
import numpy as np

## fonction de Weirstrass.
##
## n : nombre de points
## nu0 : frequence réduite de base
## nue : freq d'échantillonnage
## H : paramètre de rugosité
## choix : chaîne de caractère 
##         'd' : weierstrass déterministe
##         'a' : weierstrass aléatoire
##
## [x,t]=weir(n,nu0,nue,H,choix);

class weierstrass( object ) :
## I don't understand the mathematics behind these prameters -- must
##  investigate further!
	def __init__( self, N, H, nu0 = 1.2, nue = 1000, **kwargs ) :
		self.__N, self.__H, self.__nu0, self.__nue = N, H, nu0, nue
		self.__np_rand = None
## Reset the internal state of the generator
	def reset( self ):
		pass
## A separate procedure for lazy initialization.
	def initialize( self, numpy_random_state ) :
		self.reset( )
		self.__np_rand = numpy_random_state
## Get the approximation detail
		nmax = int( np.fix( np.log( self.__N * 0.5 ) / np.log( self.__nu0 ) ) )
## Generate time points
		self.__t = np.arange( self.__N, dtype = np.float ) / self.__nue
## Precompute values unchanged from call to call
		levels = np.arange( -nmax, nmax + 1, dtype = np.float )
		self.__scale, self.__length = self.__nu0 ** ( - levels * self.__H ), 2 * np.pi * ( self.__nu0 ** levels )
	def __call__( self ) :
## The implementation below basically constructs a fine mesh on the time axis $(t_k)_{i=0}^N\\in [0,1]$ with  
##   $$0 = t_0< t_1 < \\ldots < t_{N-1} < t_N = 1$$
##  and then evaluates the $W_\\alpha$ at the knots:
##    $$ W_H(t_i) = \sum_{k=-M}^{M} W_H^{(k)}(t_i) $$
##  with
##    $$ W_H^{(k)}( t ) = \nu_0^{-kH} \bigl( \cos(\phi_k) - \cos(2 \pi t \nu_0^k + \phi_k) \bigr) $$
##  where $M = \\bigl[\\frac{\\log \\frac{1}{2} \\nu_e}{\\log \\nu_0}\\bigr]$ and
##   $(\\phi_k)_{k=-M}^M\\sim \\mathcal{U}[0,1]$ -- random phase shift of the lay
## This array holds the values of Weierstrass function at every knot of the mesh on [0,1]
		x = np.zeros( len( self.__t ), dtype = np.float )
## Draw uniform phase shifts
		phase = 2 * np.pi * self.__np_rand.uniform( size = len( self.__scale ) )
		cos_phase = np.cos( phase )
## For each layer evaluate the trigonometric funcitons
		for k in xrange( len( phase ) ) :
			x += self.__scale[ k ] * ( cos_phase[ k ] - np.cos( self.__length[ k ] * self.__t + phase[ k ] ) )
		return self.__t, x

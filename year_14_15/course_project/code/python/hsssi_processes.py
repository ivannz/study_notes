#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from synthfbmcircul import *
from Weierstrass import *

import numpy as np

## N = power of two
def fBM( N, H = 0.6, seed = None) :
	x = synthfbmcircul( N + 1, H, sigma2 = 1.0, tmax = 1.0, seed = seed )[ 0 ]
	return x[1:]

## increasing downsample = better approximation
## Hermite process of order 2
def Rosenblatt( N, H = 0.9, seed = None, downsample = 16 ) :
## Generate fractional Gaussian Noise.
	x = synthgausscircul( N * downsample + 1, H, variance = 1.0, seed = seed )[0]
	x = np.cumsum( x**2 - 1 )[ 1: : downsample ]
## The hurst exponent for this process is H = 1 + 2 * ( Hfgn - 1 )
	return x / max( x )

## Hermite process of order 3
def Hermite3( N, H = 0.9, seed = None, downsample = 16 ) :
## Generate fractional Gaussian Noise.
	x = synthgausscircul( N * downsample + 1, H, variance = 1.0, seed = seed )[0]
	x = np.cumsum( ( x**2 - 3 ) * x )[ 1: : downsample ]
## The hurst exponent for this process is H = 1 + 3 * ( Hfgn - 1 )
	return x / max( x )

## Hermite process of order 4
def Hermite4( N, H = 0.9, seed = None, downsample = 16 ) :
## Generate fractional Gaussian Noise.
	x = synthgausscircul( N * downsample + 1, H, variance = 1.0, seed = seed )[0]
	x = np.cumsum( ( x**2 - 6 ) * x **2 + 3 )[ 1: : downsample ]
## The hurst exponent for this process is H = 1 + 4 * ( Hfgn - 1 )
	return x / max( x )

def Weier( N, H = 0.7, deterministic = False, seed = None ) :
	return synthweier( N, H = H, nu0 = 1.2, nue = 1000, seed = seed )[0]


import matplotlib.pyplot as plt
N = 2**15
t = np.arange( N, dtype = np.float ) / ( N - 1 )

a1 = fBM( N, 0.6 )
a2 = Rosenblatt( N, 0.6, downsample = 16 )
a3 = Hermite3( N, 0.9, downsample = 16 )
a4 = Hermite4( N, 0.925, downsample = 16 )
a5 = Weier( N, 0.6 )

plt.figure( figsize = (10,6) )
plt.subplot(211)
plt.plot(t, a1, "r-", linewidth = 1 )
plt.plot(t, a5, "b-", linewidth = 1 )

plt.subplot(212)
plt.plot(t, a2, "r-", linewidth = 1 )
plt.plot(t, a3, "k-", linewidth = 1 )
plt.plot(t, a4, "b-", linewidth = 1 )

plt.show( )

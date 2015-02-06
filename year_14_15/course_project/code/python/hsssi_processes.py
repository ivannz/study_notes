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

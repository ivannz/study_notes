#! /usr/bin/env python
# -*- coding: UTF-8 -*-

import re
import scipy
import numpy as np
import scipy.io as io
import scipy.stats as st
import scipy.linalg as la
# import scipy.io

mat = io.loadmat( "./Dog_1_interictal_segment_0001.mat" )

ndat = mat[ "interictal_segment_1" ]

data = ndat.item(0)[0]

l = ndat['data_length_sec'][0,0][0,0]
f = ndat['sampling_frequency'][0,0][0,0]

## Calculate running correlation matrix
def run_svd( ) :
	for j in xrange( 5, l, 5 ) :
		rho = st.spearmanr( data[:,int((j-5)*f):int((j+5)*f)], axis = 1 )
		# yield la.svd( rho[ 0 ], full_matrices = True, compute_uv = False, check_finite = False )[0]
		yield la.svdvals( rho[ 0 ], check_finite = False )[0]

np.fromiter(run_svd( ), dtype='float64')


xrange(int((j-5)*f),int((j+5)*f))

span = data.shape[1]

wnd = 1000

def run_svd( data ) :
	for t in xrange( abs(wnd), data.shape[1] ) :
		if wnd < 0 :
			rho = st.spearmanr( data[:,0:t], axis = 1 )
		else :
			rho = st.spearmanr( data[:,(t-wnd):t], axis = 1 )
		yield la.svd( rho[ 0 ], full_matrices = True, compute_uv = False, check_finite = False )[0]

## Too slow
xx = np.fromiter(run_svd(ndat.item(0)[0]), dtype='float64')


ndat.item(0)[0]

data.dot( data.transpose( ) )

## rm( list = ls( all.names = TRUE ) ) ; invisible( gc( ) )
## setwd( "~/study_notes/data_study/seizure_prediction/" )
## require( R.matlab )
## mat <- readMat( "./Dog_1_interictal_segment_0001.mat" )
## 
## ## Sampling frequency
## data <- mat$interictal.segment.1[[1]]
## 
## ## Samples in a 10 second interval
## lambda <- sapply( rev( seq( mat$interictal.segment.1[[2]], 1+5, by = -5 ) ),
## 	function( j ) {
## 		t <- seq( round( ( j - 9 ) * mat$interictal.segment.1[[3]] )
## 			, round( j * mat$interictal.segment.1[[3]] ), by = 1 )
## 		svd( cor( t( data[ ,t ] ), method = "spearman" ) )$d
## 	} )
## XX <- cor(t(mat$interictal.segment.1[[1]]), method = "spearman")
## Z<-svd(XX)python





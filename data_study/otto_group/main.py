# -*- coding: UTF-8 -*-
## Base modules
import numpy as np
import os, re

import pandas as pd

import zipfile
zf = zipfile.ZipFile( './train.csv.zip' )
df = pd.read_csv( zf.open( 'train.csv' ) )#, nrows = 1000 )

## Take a random small subsample of the train data
from sklearn.cross_validation import train_test_split

## Get the full train data and tplit it in halves
X_train_full, y_train_full = df[df.columns[1:-1]], df[df.columns[-1]]
X_train_0, X_train_1, y_train_0, y_train_1 = train_test_split( X_train_full, y_train_full, train_size = 0.50 )

#### Making an ensemble classifier
# from sklearn.metrics import log_loss
from sklearn.grid_search import GridSearchCV
def scores_to_df( grid ) :
	res = pd.DataFrame( par[ 0 ] for par in grid )
	res[ 'metric' ] = pd.DataFrame( par[ 1 ] for par in grid )
	return res

## The strategy is to use as many classifiers as possible in hope that each one could 
##  capture some hidden structural aspects of the data so that together the classifiers 
##  would have better performance.

ensemble_clf = list( )

## 1. LogisticRegression
from sklearn.linear_model import LogisticRegression
logreg_grid = GridSearchCV( LogisticRegression( multi_class = "ovr" ), param_grid = {
		"C" : np.logspace( -1, 2, num = 4 ),
	}, cv = 10, n_jobs = -1, verbose = 50, scoring = "log_loss" ).fit( X_train_0, y_train_0 )

ensemble_clf.append( ( "Logistic", logreg_grid.best_estimator_ ) )
scores_to_df( logreg_grid.grid_scores_ )

## 2.Random Forest
from sklearn.ensemble import RandomForestClassifier
rf_grid = GridSearchCV( RandomForestClassifier( n_estimators = 256 ), param_grid = {
## max_depth -- the maximum allowed number of levels in the decision tree.
		"max_depth" : [ 1, 3, 5, 7, 10, ],
	}, cv = 10, scoring = 'log_loss', verbose = 10 ).fit( X_train_0, y_train_0 )

ensemble_clf.append( ( "Forest", rf_grid.best_estimator_ ) )
scores_to_df( rf_grid.grid_scores_ )

## 3. k - nearest neighbour
from sklearn.neighbors import KNeighborsClassifier
knn_clf = KNeighborsClassifier( n_neighbors = 2 ).fit( X_train_0, y_train_0 )
## Survey the data landscape with the Nearest Neighbours Classifiers
knn_clf = [
	( "knn-%d" % ( n_neighbors, ), KNeighborsClassifier( n_neighbors = n_neighbors ).fit( X_train_0, y_train_0 ) )
		for n_neighbors in [ 2, 8, 32, 128, 512, ] ]

ensemble_clf.extend( knn_clf )

## 4. XGboost
from xgboost import XGBClassifier
xgb_clf = XGBClassifier( )



## show log loss
from sklearn.preprocessing import LabelBinarizer
lb = LabelBinarizer( ).fit( df[ df.columns[ -1 ] ] )
y_test_ovr = lb.transform( y_train_1 )

for name, clf in ensemble_clf :
	theta = clf.predict_proba( X_train_1 )
	A = np.tensordot( y_test_ovr, np.log( np.clip( theta, 1e-15, 1-1e-15) ), ( 0, 0 ) )
	logloss = -np.sum( np.diag( A ), dtype = np.float ) / y_test_ovr.shape[ 0 ]
	print "%s: logLoss %.5f" % ( name, logloss, )

## Make a cube of predicted probabilities NxKxC
proba_cube = np.zeros( ( y_train_1.shape[0], len( lb.classes_ ), len( ensemble_clf ), ), dtype = np.float )
for k, ( name, clf ) in enumerate( ensemble_clf[:2] ) :
	proba_cube[:,:,k] = clf.predict_proba( X_train_1 )
	pass

print proba_cube.sum( axis = 1 )

## Combine the classifiers by a uniform weight
weight = np.asarray( [ .5, .5 ], dtype = np.float )
theta = np.tensordot( weight, proba_cube, (0, 2) )
print -np.sum( np.diag( A ), dtype = np.float ) / y_test_ovr.shape[ 0 ]

####################################################################################################

# from sklearn.metrics import log_loss
## ## Encode the class assignments as a 1-of-all vector
## from sklearn.preprocessing import LabelBinarizer
## lb = LabelBinarizer( ).fit( df[ df.columns[ -1 ] ] )
## ## Compute the logloss
## e_hat = lb.transform( y_hat )
## A = np.tensordot( e_hat, np.log( proba ), ( 0, 0 ) )
## logloss = -np.sum( np.diag( A ) / e_hat.shape[ 0 ] )

## Predict the class labels
proba = lreg.predict_proba( X_test )

## Compute the logloss
e_test = lb.transform( y_test )

## Logloss = -\frac{1}{N} \sum_{i=1}^N \sum_{k=1}^K t_{ik} \log \hat{p}_{ik}
## 		   = - \sum_{k=1}^K  \frac{1}{N} A_{kk}, where
## 	A_{jk} =\sum_{i=1}^N t_{ij} \log \hat{p}_{ik}
A = np.tensordot( e_test, np.log( proba ), ( 0, 0 ) )
logloss = -np.sum( np.diag( A ) / e_test.shape[ 0 ] )

####################################################################################################




## Encode the class assignments as a 1-of-all vector
from sklearn.preprocessing import LabelBinarizer
lb = LabelBinarizer( ).fit( df[ df.columns[ -1 ] ] )

## Train a set of logixtic regressions with 1-vs-ALL classification
from sklearn.linear_model import LogisticRegression
lreg = LogisticRegression( multi_class = 'ovr' ).fit( X_train_full, y_train_full )

## Load the test sample
zf = zipfile.ZipFile( './test.csv.zip' )
df_test = pd.read_csv( zf.open( 'test.csv' ) ) #, nrows = 100 )

## Get the independent variables
X_test = df_test[df.columns[1:-1]]

proba = lreg.predict_proba( X_test )
y_hat = lreg.predict( X_test )

## Compute the logloss
e_hat = lb.transform( y_hat )
A = np.tensordot( e_hat, np.log( proba ), ( 0, 0 ) )
logloss = -np.sum( np.diag( A ) / e_hat.shape[ 0 ] )


result = pd.DataFrame( { df.columns[ 0 ] : df_test[ df.columns[ 0 ] ] } )
result[ lb.classes_ ] = pd.DataFrame( { k : e_hat[:,j] for j, k in enumerate( lb.classes_, 0 ) } )
result.to_csv( "./sampleSubmission.csv", index = False )


####################################################################################################

from sklearn.ensemble import RandomForestClassifier

# from sklearn.ensemble import AdaBoostClassifier

from sklearn.grid_search import GridSearchCV

rf_grid = GridSearchCV( RandomForestClassifier( n_estimators = 256 ), cv = 10,
		param_grid = { "max_depth": [ 3, 5, 12, 25, ], }, verbose = 10, n_jobs = -1 ).fit( X_train, y_train )

clf = rf_grid.best_estimator_.fit( X_train, y_train )

y_hat = clf.predict( X_test )

##################
zf = zipfile.ZipFile( './test.csv.zip' )
df_test = pd.read_csv( zf.open( 'test.csv' ) ) #, nrows = 100 )

y_test_hat = clf.predict( df_test[df.columns[1:-1]] )

result = pd.DataFrame( { df.columns[ 0 ] : df_test[ df.columns[ 0 ] ] } )
class_labels = np.unique( df[ "target" ].values )
result[ class_labels ] = pd.DataFrame( { k : (y_test_hat == k)*1 for k in class_labels } )

result.to_csv( "./out.csv", index = False )

###################
## EM algorithm ##

import scipy as sp

theta_init = np.random.gamma( 1.0, size = ( X_train.shape[ 1 ], 9 ) )
theta_init /= np.sum( theta_init, axis = 0 ).reshape( ( 1, -1 ) )

## 
sp.special.gammaln( .5 )

####################################################################################################
from sklearn.preprocessing import LabelBinarizer
lb = LabelBinarizer( ).fit( df[df.columns[-1]] )

## Convert to 1-of-all encoding
t_train = lb.transform( y_train )
t_test = lb.transform( y_test )

class_freq = np.tensordot( X_train, t_train, ( 0, 0 ) )
pi = np.sum( t_train, axis = 0, dtype = np.float ).reshape( ( 1, -1 ) ) / t_train.shape[ 0 ]

################################# (Categorical) Multinomial model #################################
theta = class_freq / np.sum( class_freq, axis = 0, dtype = np.float ).reshape( ( 1, -1 ) )

ll_mult = np.dot( X_test, np.log( theta ) ) + np.log( pi )
t_hat_mult = np.argmax( ll_mult, axis = 1 )

################################# (Categorical) Multinomial model #################################
feature_rate = class_freq / np.sum( t_train, axis = 0, dtype = np.float ).reshape( ( 1, -1 ) )

ll_pois = np.dot( X_test, np.log( feature_rate ) ) - np.sum( feature_rate, axis = 0 ) + np.log( pi )
t_hat_pois = np.argmax( ll_pois, axis = 1 )

cl_hat_pois = lb.classes_[t_hat_pois]

np.sum( cl_hat_pois == y_test )

####################################################################################################
np.unique( t_hat_mult, return_counts = True )
np.unique( t_hat_pois, return_counts = True )

####

from sklearn.neighbors import KNeighborsClassifier
from sklearn.grid_search import GridSearchCV

knn_grid = GridSearchCV( KNeighborsClassifier( ), cv = 10, n_jobs = -1, verbose = 10,
	param_grid = { "n_neighbors" : [ 1, 3, 5, 12, 25, 60 ] } ).fit( X_train, y_train )

knn_clf = knn_grid.best_estimator_.fit( X_train, y_train )

y_hat = knn_clf.predict( X_test )

## Make the confusion matrix
from sklearn.preprocessing import LabelEncoder
le = LabelEncoder( ).fit( df[df.columns[-1]].values )

i_test, i_hat = le.transform( y_test ), le.transform( y_hat )

confusion = np.zeros( 2 * [ len( le.classes_ ) ], dtype = np.int )
for i in range( confusion.shape[0] ) :
	j, f = np.unique( i_test[ i_hat == i ], return_counts = True )
	confusion[ i, j ] = f

print confusion


from sklearn.pipeline import Pipeline, FeatureUnion

## Survey the data landscape with the Nearest Neighbours Classifiers
knn_clf = [
	( "knn-%d" % ( n_neighbors, ), KNeighborsClassifier( n_neighbors = n_neighbors ).fit( X_train, y_train ) )
		for n_neighbors in [ 2, 4, 8, 16, 32, 64, 128, 256 ] ]

nn_features = pd.DataFrame( { name : knn.predict( X_train ) for name, knn in knn_clf } )

X_train[nn_features.columns] = nn_features


pipeline = Pipeline( [
	( "features", combined_features ),
## Use rbf-kernel svm
	# ( "svm", SVC( kernel = "linear" ) ),
	( "svm", SVC( kernel = "rbf" ) ),
	# ( "forest", RandomForestClassifier( ) ),
] )



le = LabelEncoder( ).fit( df[df.columns[-1]].values )

i_train, i_hat = le.transform( y_train ), le.transform( nn_features.values[:,8] )

confusion = np.zeros( 2 * [ len( le.classes_ ) ], dtype = np.int )
for i in range( confusion.shape[0] ) :
	j, f = np.unique( i_train[ i_hat == i ], return_counts = True )
	confusion[ i, j ] = f

print confusion






import numpy as np
import scipy.sparse as sp

A = sp.csr_matrix( [
	[0,1,0,0,1,0],
	[0,0,1,0,0,1],
	[0,0,0,1,1,0],
	[0,0,0,0,1,0],
	[0,1,0,0,0,0],
	[0,0,0,0,0,0], ], shape = ( 6, 6 ), dtype = np.float )

## Out-degree
out_deg = np.asarray( A.sum( axis = 1 ).getA1( ), dtype = np.float )
dangling = np.where( out_deg == 0 )[ 0 ]
out_deg[ dangling ] = 1.0

beta = 0.85
E = np.full( A.shape[ 0 ], 1.0 / A.shape[ 0 ], dtype = np.float )

## Since matrix is a linear operator and the eigenvalues we seek is one, the requirement
##  that the scores vector sum to one is automatically stisfied once it has been imposed.
x_0 = E.copy( )
x_0 = beta * ( x_0 / out_deg ) * A + beta * np.sum( x_0[ dangling ] ) * E + ( 1 - beta ) * np.sum( x_0[ 0 ] ) * E

## D = \text{diag}\bigl( \delta^+_v + 1_{\delta^+_v=0} \bigr)_{v\in V}\,,
## d = \bigl( 1_{\delta^+_v=0} \bigr)_{v\in V}\,, the indicator vector of dangling vertices
##  Q = D^{-1} A + d ( 1'1 )^{-1} 1'\,,
## The matrix with teleportation option:
##  M = \beta Q + ( 1 - \beta ) 1 ( 1'1 )^{-1} 1'\,.
## For a personalized pagerank for w\in V use e_w = ( 1_{v=w} )_{v\in V}\,,
##  ... + ( 1 - \beta ) 1 (1' e_w)^{-1} e_w'\,.

## Y = X M = \beta X Q + (1-\beta) X 1 (1'1)^{-1} 1' 
## 		   = \beta X D^{-1} A + \beta X d ( 1'1 )^{-1} 1' + (1-\beta) X 1 (1'1)^{-1} 1' 


def sparse_pagerank( A, beta = 0.85, one = None, niter = 1000, rel_eps = 1e-6 ) :
## Initialize the iterations
	one = one if one is not None else np.ones( ( 1, A.shape[ 0 ] ), dtype = np.float )
	one = sp.csr_matrix( one / one.sum( axis = 1 ) )
## Get the out-degree
	out = np.asarray( A.sum( axis = 1 ).getA1( ), dtype = np.float )
## Obtain the mask of dangling vertices
	dangling = np.where( out == 0.0 )[ 0 ]
## Correct the out-degree for sink nodes
	out[ dangling ] = 1.0
## Just one iteration: all dangling nodes add to the importance of all vertices.
	pi = np.full( ( one.shape[0], A.shape[0] ), 1.0 / A.shape[ 0 ], dtype = np.float )
## If there are no dangling vertices then use simple iterations
	kiter, status = 0, -1
## Make a stochastic matrix
	P = sp.diags( 1.0 / out, 0, dtype = np.float ).dot( A ).tocsc( )
	while kiter < niter :
## make a copy of hte current ranking estimates
		pi_last = pi.copy( )
## Use sparse inplace operations for speed. Firstt the random walk part
		pi *= beta ; pi *= P
## Now the teleportaiton ...
		pi += ( 1 - beta ) * one
##  ... and dangling vertices part
		if len( dangling ) > 0 :
			pi += beta * one.multiply( np.sum( pi_last[ :, dangling ], axis = 1 ).reshape( ( -1, 1 ) ) )
## Normalize
		pi /= np.sum( pi, axis = 1 )
		if np.sum( np.abs( pi - pi_last ) ) <= one.shape[0] * rel_eps * np.sum( np.abs( pi_last ) ) :
			status = 0
			break
## Next iteration
		kiter += 1
		if kiter % 10 == 0 :
			print kiter
	return pi, status, kiter


pi1, s, k = sparse_pagerank( A, one = None )

one = sp.eye( A.shape[ 0 ] ).tocsr()[:100]
pi2, s, k = sparse_pagerank( A, one = one )

one = sp.eye( A.shape[ 0 ] ).tocsr()[:100]
pi3, s, k = sparse_pagerank( A, one = one )

x0 = np.full( A.shape, 1.0 / A.shape[ 0 ], dtype = np.float )
for i in range( 550 ):
	x0 = beta * ( x0 / out ) * A + ( beta * np.sum( x0[ :, dangling ], axis = 1 ).reshape( ( -1, 1 ) ) + ( 1 - beta ) ) * one

## One iteration of the power iterations method
	x1 = beta * ( x0 / out ) * A + ( beta * np.sum( x0[ dangling ], axis =  ) +  ( 1 - beta ) ) * one


	x1 = beta * ( x0 / out ) * A + ( beta * np.dot( x0.T, d ) + ( 1 - beta ) * np.dot( x0.T, one ) ) * one.T
	pass


import networkx as nx
G = nx.from_scipy_sparse_matrix( A, create_using = nx.DiGraph( ) )
nx.pagerank_scipy( G )


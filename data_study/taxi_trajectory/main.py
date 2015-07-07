import json
import zipfile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time as tm
from ast import literal_eval as from_WGS84

# reading training data
tick = tm.time( )
zf = zipfile.ZipFile( './data/train.csv.zip' )
df = pd.read_csv( zf.open( 'train.csv' ), # nrows = 50000,
	converters = { 'POLYLINE' : lambda x : np.asarray( json.loads( x ) ).reshape( -1, 2 ) } )
tock = tm.time( )
print "Taxi trip data loaded in %.3f sec." % ( tock - tick, )
latlong = np.concatenate( df[ 'POLYLINE' ].values, axis = 0 )

# df_old = pd.read_csv( zf.open( 'train.csv' ), converters = { 'POLYLINE' : lambda x : json.loads( x )[ -1 : ] } )
# latlong = np.array( [[ p[ 0 ][ 1 ], p[ 0 ][ 0 ] ] for p in df[ 'POLYLINE' ] if len( p ) > 0 ] )

# cut off long distance trips
lat_low, lat_hgh = np.percentile( latlong[ :, 0 ], [ 2, 98 ] ) 
lon_low, lon_hgh = np.percentile( latlong[ :, 1 ], [ 2, 98 ] )

# create image
bins = 513
lat_bins = np.linspace( lat_low, lat_hgh, bins )
lon_bins = np.linspace( lon_low, lon_hgh, bins )
H2, _, _ = np.histogram2d( latlong[ :, 0 ], latlong[ :, 1 ], bins = ( lat_bins, lon_bins ) )

img = np.log( H2[ : : -1, : ] + 1 )

plt.figure( )
ax = plt.subplot( 1, 1, 1 )
plt.imshow( img )
plt.axis( 'off' )
plt.title( 'Taxi trip end points' )
plt.show( )

plt.savefig( "taxi_trip_end_points.png" )



from sklearn.cluster import DBSCAN

from math import radians, cos, sin, atan, sqrt

## The callable should take two arrays from X as input and return a value indicating
##  the distance between them.
def haversine_( X1, X2 ) :
	print X1.shape, X2.shape
	(lat1, lon1), (lat2, lon2) = X1, X2
	a = sin( ( lat2 - lat1 ) / 2 )**2 + cos( lat1 ) * cos( lat2 ) * sin( ( lon2 - lon1 ) / 2 )**2
## 2 r = 2 * 6378137
	return 12756274 * atan( sqrt( a / ( 1 - a ) ) )

# Compute DBSCAN
db = DBSCAN( eps = 15, min_samples = 10, metric = haversine_ ).fit( latlong[:1000] )

core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
core_samples_mask[ db.core_sample_indices_ ] = True
labels = db.labels_

# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

# Black removed and is used for noise instead.
unique_labels = set(labels)
colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
for k, col in zip( unique_labels, colors ) :
    if k == -1:
        # Black used for noise.
        col = 'k'

    class_member_mask = (labels == k)

    xy = latlong[ class_member_mask & core_samples_mask ]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=14)

    xy = latlong[ class_member_mask & ~core_samples_mask ]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=6)

plt.title('Estimated number of clusters: %d' % n_clusters_)
plt.show()

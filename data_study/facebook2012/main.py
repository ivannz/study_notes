#! /usr/bin/env python

from collections import deque
import utilities as uio

import base as base
import numpy as np

## Load the graph
adj = base.__csr_from_csv( "train.gz", compression = "gzip" )

## Path shortening 
pi = base.__sparse_bfs( adj, [1,], 100, 2 )
print np.where( ( pi > 1 ) & np.isfinite( pi ) )[0].__len__()

## target nodes
target = pd.read_csv( "test.csv" ).values[ :, 0 ]

def common_neighbours( adj ) :
	pass

####################################################################################################
################################ Old code from the provided example ################################
####################################################################################################

def breadth_first_search( graph, node, num_nodes ) :
	"""
	Does a breadth-first search of the graph starting at the node.
	Returns the first num_nodes nodes (excluding direct neighbors)
	"""
	neighbors = set( graph[ node ] )
	looked_at = set( graph[ node ] )
	looked_at.add( node )
	visited = [ ]
	queue = deque( graph[ node ] )

	while queue and len( visited ) < num_nodes :
		next_node = queue.popleft( )
		if next_node not in neighbors :
			visited.append( next_node )
		queue.extend( n for n in graph[ next_node ] if n not in looked_at )
		looked_at.update( graph[ next_node ] )

	return visited

def bfs_benchmark(train_file, test_file, submission_file, num_predictions):
	"""
	Runs the breadth-first search benchmark.
	"""
	graph = utilities.read_graph(train_file)
	test_nodes = utilities.read_nodes_list(test_file)
	test_predictions = [breadth_first_search(graph, node, num_predictions)
						for node in test_nodes]
	utilities.write_submission_file(submission_file, 
									test_nodes, 
									test_predictions)

if __name__=="__main__":
	bfs_benchmark("train.csv.gz",
				  "test.csv",
				  "../Submissions/bfs_benchmark_updated.csv",
				  10)
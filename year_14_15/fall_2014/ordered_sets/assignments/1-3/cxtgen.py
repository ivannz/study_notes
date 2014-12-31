#! /usr/bin/env python
# -*- coding: UTF-8 -*-
import random
import sys

def usage() :
	print "%s |G| |M| context" % (__file__)
	return -1

if __name__ == '__main__' :
	if len( sys.argv ) < 3:
		exit( usage() )

	G = int( sys.argv[ 1 ] )
	M = int( sys.argv[ 2 ] )
	filename_out = sys.argv[ 3 ]

	lines = list( )
	lines.append( "B\ntest\n%d\n%d\n"%(G,M) )
	lines.append( "\n".join( str( x ) for x in xrange( G ) ) + "\n" )
	lines.append( "\n".join( str( x ) for x in xrange( M ) ) + "\n" )
	for j in xrange( G ) :
		x = random.getrandbits( M )
		lines.append( "".join( ".X"[1 & ( x >> i )] for i in xrange( M ) ) + "\n" )

	with file( filename_out, "w+" ) as f:
		f.writelines( lines )
	f.close()


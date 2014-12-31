#! /usr/bin/env python
# -*- coding: UTF-8 -*-
####################################################################################################
####################################################################################################
####################################################################################################
##### USAGE: python cbo.py context output intent/extent-flag
##### Formal concept analysis begins at line 106. Look for ##FCA_HERE## marker.
#####  WARNING: THE CODE IS HEAVILY COMMENTED!
__date__ 			= "2014-10-30"
__author__ 			= "Nazarov Ivan"
__email__ 			= "innazarov@edu.hse.ru"
__status__ 			= "Barely useable"
__version__ 		= "0.8"
__dscription__ 		= "Python realisation of CbO with the possibility of fc lattice recovery. Assigment #03."

"""Sorry, no Russian while coding!"""

####################################################################################################
####################################################################################################
####################################################################################################
import re
import sys
###### intbitset is used for fast hardware optimized set operations
###      http://intbitset.readthedocs.org/en/latest/#
import intbitset as bs

## This tiny function calculates a unique identifier of a subset of integers.
def label( S ) :
	i = 0
	for x in S : i += 1 << x
	return i

## a Very basic class to store forward cascading connectivity
class successor_node(object):
	def __init__(self, tag, pair):
		self.tag = tag
		self.pair = pair
		self.__prev = None
		self.__next = set()
	def alpha(self):
		while self.__prev is not None:
			self = self.__prev
		return self
	def generator(self) :
		return self.__prev
	def successor(self) :
		return frozenset( self.__next )
	def link(this, that) :
		if isinstance( that, successor_node ) :
			that.detach()
			this.__next.add( that )
			that.__prev = this
		return that
	def detach(self) :
## Disconnect from the  parent
		prev = self.__prev
		if isinstance( prev, successor_node ) :
			prev.__next.remove( self )
		return prev
	def sever(self) :
## Sever forward links
		for next in self.__next :
			if isinstance( next, successor_node ) and next.__prev is self :
				next.__prev = self.__prev
		return self
	def unlink(self) :
		self.detach()
		self.sever()
## Clean up
		self.__next.clear()
		self.__prev = None
		return self
	# def __repr__(self) :
		# return "-".join([str(e) for e in self.pair])

## a Very basic class to store lattice connectivity
class lattice_node(object):
	def __init__(self, tag, pair):
		self.tag = tag
		self.pair = pair
		self.sup = set()
		self.inf = set()
	def link_sup(this, that) :
		if isinstance( that, lattice_node ) :
			this.sup.add( that )
			that.inf.add( this )
	def link_inf(this, that) :
		if isinstance( that, lattice_node ) :
			this.inf.add( that )
			that.sup.add( this )
	def unlink_inf(self):
		for i in self.inf : i.sup.remove(self)
		self.inf.clear()
	def unlink_sup(self):
		for s in self.sup : i.inf.remove(self)
		self.sup.clear()
	def get_inf(self) :
		while len( self.inf ) > 0:
			self = next( iter( self.inf ) )
		return self
	def get_sup(self) :
		while len( self.sup ) > 0:
			self = next( iter( self.sup ) )
		return self
	def unlink(self) :
		self.unlink_inf()
		self.unlink_sup()
	# def __repr__(self) :
		# return "-".join([str(e) for e in self.pair])

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
""" ##FCA_HERE## Formal Context analysis"""
def skip_empty( f ):
# Skip an empty line if there is one
	pos = f.tell()
	if f.readline().rstrip():
		f.seek(pos)

def burmeister(filename) :
	with file( filename, "r" ) as f:
# File format marker
		txt = f.readline().strip();
		if txt is not "B":
			raise Exception("Bad file format")
# The name of the context
		name = f.readline().strip();
# Read the volume of the extent and the intent
		G = int(f.readline().strip());
		M = int(f.readline().strip());
		skip_empty( f )
# Read the labels first the extent, last the intent
		extent = [ f.readline().strip() for x in range( G ) ]
		intent = [ f.readline().strip() for x in range( M ) ]
		skip_empty( f )
# Load and process the object-attribute relationship,
#   expect that it is stored extent-wise
		ctx0 = []
		for g in range( G ):
			line = f.readline().rstrip()
			ctx0.append(bs.intbitset( [
				m for m in range( M ) if line[ m ] in "X" ], M ) )
# Transpose the context
		ctxt = [ bs.intbitset( [ g for g in range( G )
			if m in ctx0[ g ] ], G ) for m in range( M ) ]
		return Context( extent, intent, ( ctx0, ctxt ), name )

class Context(object):
	def __init__(self, extent, intent, ctx, name):
		super(Context, self).__init__()
		self.intent = intent
		self.extent = extent
		self.name = name
## Prepare the context
		self.__ctx = ctx
		self.__extent = bs.intbitset( range( len( self.extent ) ) )
		self.__intent = bs.intbitset( range( len( self.intent ) ) )
	def G(self, S):
## return the labels of the extent
		return [ self.extent[ x ] for x in bs.intbitset( S ) ]
	def M(self, S):
## return the labels of the intent
		return [ self.intent[ x ] for x in bs.intbitset( S ) ]
	def GM(self, S, dual):
## Represent the set of volumes 
		if dual :
			return [ ( "|".join( self.G( c[ 1 ] ) ), "|".join( self.M( c[ 0 ] ) ) ) for c in S ]
		else :
			return [ ( "|".join( self.G( c[ 0 ] ) ), "|".join( self.M( c[ 1 ] ) ) ) for c in S ]
################################################################################
## The context (G,M,R) defines a Galois connection between the power sets of G
##  and M. u:2^G->2^M and d:2^G->2^M possess the following properties:
##  for all A,X ≤ G and B,Y ≤ M the following are true:
##    GALOIS: 	A ≤ dB <=> B ≤ uA
##    AMON: 	A ≤ X => uX ≤ uA ; B ≤ Y => dY ≤ dB
##    UNION: 	u(A|X) = uA & uX ; d(B|Y) = dB & dY
##  Different composition of "u" and "d" define closure operators on G and M:
##    CLOS: 	a) A ≤ duA; b) A≤C => duA ≤ duC; c) duduA = duA
##    			a) B ≤ udB; b) B≤D => udB ≤ udD; c) ududB = udB
	def u(self, S) :
## The "up" function maps a subset of G to a subset of M through the context
##  relation R
		res = self.__intent
		return set(res.intersection( *[ self.__ctx[ 0 ][ x ] for x in bs.intbitset( S ) ] ))
	def d(self, S) :
## THe "down" map returns a subset of G which corresponds to the given subset of M
		res = self.__extent
		return set(res.intersection( *[ self.__ctx[ 1 ][ x ] for x in bs.intbitset( S ) ] ))
	# def reduce(self):
	# 	pass
	def __CbO(self, P, G, u, d):
## This function initiates the CbO for the Galois connection (u,d)
##  In future one should abstract the Galois connection into a separate class.
## For any A ≤ G: uA = uduA; since A ≤ duA, AMON implies uduA ≤ uA whereas
##  by GALOIS from duA ≤ duA follows uA ≤ uduA. Therefore the following
##  invocation to __CbO_step as it initiates the Close-by-one algorithm form
##  a correct dual pair (duP, uP), where P is considered to be the seed set.
		suxx = self.__CbO_step( d(u(P)), u(P), G, u, d, 0, {}, 0 )
## Construct a lattice
		ltce = self.__lattice( d(u(P)), u(P), G, u, d, {} )
		return (suxx, ltce)
################################################################################
################################################################################
################################################################################
	def __CbO_step(ctx, A, B, G, u, d, g, suc, lvl):
################################################################################
## Based upon multiple sources:
##  Kuznetsov, S.; Obiedkov, S (2002).
##     Comparing Performance of Algorithms for Generating Concept Lattices.
##     Journal of Experimental & Theoretical Artificial Intelligence, 14.
##     doi:10.1.1.4.1110
##  Krajca P., Outrata J., Vychodil V.
##     Advances in algorithms based on CbO
##     Presentation at Palacky University, Olomouc, Czech Republic
##         Concept Lattices and Their Applications, Sevilla, 2010
################################################################################
## This is a recursive step of a CbO algorithm for the connection (u,d)
##  Expect (A,B) to be a dual pair with respect to the supplied Galois
##   connection: B = uA, A = dB, where u:Ω_A -> Ω_B is the Galois map
##   to obtain the reflection of A in the universe of B and d:Ω_B -> Ω_A
##   reflects B into A's universe. G is the size of A's universe (|Ω_A|)
##DEBUG##		assert( B == u(A) ) ; assert( A == d(B) )
## This function is called only when a new concept has been canonically generated
## The lattice is indexed by the spawning set
		suc[ label( A ) ] = successor_node( label( A ), ( A, B ) )
## Examine elements in lexicographic order unless there is nothing more to do
		for f in range( g, G ) : # for f in set( range( g, G ) ) - A :
## There is no point in trying to expand by a member element, since we will generate
##  ourselves again
			if f in A : continue
## Close-by-One: generate a concept by adding a new member to it:
##   it is true that u(A | {f}) = uA & u{f} = B & u{f} (UNION)
			Z = B & u( { f } ) ; Y = d( Z )
##DEBUG##			print " "*lvl, "%d»» %s %d]%s)" % ( lvl, "".join([str(e) for e in A]), f, "".join([str(g) for g in (Y-(A|{f}))]) )
## Perform the canonicity test: if all elements of Y \ A ( the newly generated
##  volume) are greater of equal to f then the generation (Y,Z) from A is
##  canonical. The test min( h | h in du(A | {f})-A ) ≥ f is actually equivalent
##  to checking if (du(A | {f})-A) & F is ø, where F={0..f-1}. In turn, since
##  C ≤ duC (CLOS) this test is equivalent to checking if A & F = du(A | {f}) & F
			F = set( range( f ) )
			if F & Y == F & A :
## Every canonically generated pair (Y,uY) is guaranteed to have never been
##  generated by anyone before, because pairs are generated in lexicographic order.
##  Generate from the newly detected canonical dual pair (Y,uY)
				ctx.__CbO_step( Y, Z, G, u, d, f + 1, suc, lvl + 1 )
## Add the newly generated dual pair to the list of descendants of the current pair
				suc[ label( A ) ].link( suc[ label( Y ) ] )
		return suc
################################################################################
################################################################################
################################################################################
	def __lattice(ctx, A, B, G, u, d, L) :
################################################################################
## A simplified and tested version of [http://ami.hse.ru/data/165/315/1234/book.pdf] p. 117
################################################################################
## Create a new lattice node
		L[ label( A ) ] = lattice_node( label( A ), ( A, B ) )
## Collect all dual pairs that are generated from the current one by adding
##  a single new element 
		queue = dict()
		for j in range( G ) :
## Again, as in CbO above, there is no use in adding a part of the existing
##  dual pair
			if j in A : continue
## Like in CbO generate the next pair by adding a new element and closing with
##  the Galois connection (u,d)
			Y = B & u( { j } ) ; X = d( Y )
##DEBUG##		print "»» %s %d]%s)" % ( "".join([str(e) for e in A]), j, "".join( [str(g) for g in ( Z - ( A | { j } )  ) ] ) )
## Add the new heir to the list
			queue[ label( X ) ] = ( X, Y )
## Weed out nested heirs unfortunately due to the fact that set inclusion is 
##  inly a partial order on subsets, this has cubic complexity. This constructs
##  the next layer of nodes, which are related to the current via a covering
##  relation based on lattice partial order
		next = dict()
		for M in queue.values( ) :
## Scan the siblings for the smallest one with respect to the current M
			for Z in queue.values( ) :
				if Z[0] < M[0] : M = Z
## Add it to the call list.
			next[ label( M[0] ) ] = M
## Recursive call:
		for X,Y in next.values( ) :
## If the successor has already been visited, there is no need to traverse
##  its sub-lattice
			if label( X ) not in L :
				ctx.__lattice( X, Y, G, u, d, L )
## Connect this successor to the current node
			L[ label( A ) ].link_sup( L[ label( X ) ] )
		return L
	def CbOG(self):
## Generate by closing extents # (A,B) if not dual else (B,A)
		result = self.__CbO( set(), len( self.extent ), self.u, self.d )
		return ( result, False )
	def CbOM(self):
## Generate by closing intents # (A,B) if not dual else (B,A)
		result = self.__CbO( set(), len( self.intent ), self.d, self.u )
		return ( result, True )
	def save_lattice( self, L, filename, dual ):
		if not L :
			return
		lines = list()
## Pick a random key from the lattice dictionary and from it find the root pair
		R = L[ L.keys()[ 0 ] ].get_inf( )
## Get the number of pairs and edges
		n = len( L )
		e = sum( [ len(x.sup) for x in L.values() ] )
		lines.append( "%d\n%d\n" %( n, e ) )
## Identify pairs with consecutive natural numbers
		np = L.keys() ; pn = dict()
		for x in range( n ) :
			pn[ np[ x ] ] = x
			A, B = self.GM( [ L[ np[ x ] ].pair ], dual )[ 0 ]
			lines.append( "%d '%s' '%s'\n" % ( x, A, B ) )
## Recover adjacency lists
		for k, x in L.items() :
			if x.sup :
				lines.append( "%d %s\n" %( pn[ k ], " ".join( [ str( pn[ z.tag ] ) for z in x.sup ] ) ) )
		with file( filename, "w+" ) as f:
			f.writelines( lines )

def usage() :
	print "%s <path-to-context> <output-filename> <intent/extent>" % (__file__)
	return -1

if __name__ == '__main__' :
	if len( sys.argv ) < 3:
		exit( usage() )

	filename_in = sys.argv[1]
	filename_out = sys.argv[2]
	intent_driven = True 
	if len(sys.argv) == 4 :
		intent_driven = True if sys.argv[3].lower() == 'y' else False

	ctx = burmeister(filename_in);

	if intent_driven :
## Intent driven CbO
		(Tree, Lattice), Duality_Flag = ctx.CbOM()
	else :
## Extent driven CbO
		(Tree, Lattice), Duality_Flag = ctx.CbOG()

## Save the lattice
	ctx.save_lattice( Lattice, filename_out, Duality_Flag )

#! /usr/bin/env python
# -*- coding: UTF-8 -*-
####################################################################################################
####################################################################################################
####################################################################################################
##### USAGE: python cbo.py context output intent/extent-flag
##### Formal concept analysis begins at line 106. Look for ##FCA_HERE## marker.
#####  WARNING: THE CODE IS HEAVILY COMMENTED!
__date__ 			= "2014-12-27"
__author__ 			= "Nazarov Ivan"
__email__ 			= "innazarov@edu.hse.ru"
__status__ 			= "Useable"
__version__ 		= "1.2"
__dscription__ 		= "Python realisation of CbO on pattern structures with the possibility of pattern lattice recovery. Assigment #01 module #02."

"""Sorry, no Russian while coding!"""

####################################################################################################
####################################################################################################
####################################################################################################
import sys

## This tiny function calculates a unique identifier of a subset of integers.
def label( S ) :
	i = 0
	for x in S : i += 1 << x
	return i

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
	def __repr__( self ) :
		return "-".join( [ str( e ) for e in self.pair ] )

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
## Define structures for semilattice with join-semantics
## Algebraic: operarion d1 | d2 is guaranteed to be idempotent,
##  associative and commutative by construction.

## Implements descriptions which are subsets of some finite set 
class S(object) :
	def __init__(self, volume = []) :
		super(S, self).__init__()
		self.volume = frozenset( volume )
	def __or__( i, j ) :
### The join operation on two descriptions must produce a more
###  general description. The notation goes against the one
###  usd in lectures
		return S( i.volume | j.volume )
	def __le__( i, j ) :
### Must be equivalent to j == i | j
		return i.volume <= j.volume
	def __eq__( i, j ) :
		return i.volume == j.volume
	def __ne__( i, j ) :
		return i.volume != j.volume
	def __repr__( self ) :
		if not self.volume : return "ø"
		return "{%s}" % ",".join([str(e) for e in self.volume])

## Implements interval descriptions suitable for continuous numeric features
class I(object) :
	def __init__(self, lim = []) :
		super(I, self).__init__()
		self.lim = tuple( lim )
	def __or__( i, j ) :
### Lazy interval envelope: the smallest interval covering both
		return I( i.lim + j.lim )
	def __le__( i, j ) :
### Must be equivalent to j == ( i | j ) (join-semantics)
		if not i.lim : return True
		if not j.lim : return False
		i.__normalize( ) ; j.__normalize( )
		return j.lim[0] <= i.lim[0] and i.lim[1] <= j.lim[1]
	def __eq__( i, j ) :
		j.__normalize( ) ; i.__normalize( )
		return i.lim == j.lim
	def __ne__( i, j ) :
		return not (i == j)
	def __repr__( self ) :
		if not self.lim : return "ø"
		self.__normalize( )
		return "[%g,%g]" % self.lim
	def __normalize( self ) :
		if self.lim :
			self.lim = tuple( [ min( self.lim ), max( self.lim ) ] )

## Implement a semilattice product
## I am very sorry Sergei, but join semantics in description
##  semilattices is much more intuitive for me.
class D(object) :
	def __init__( self, volume = [] ) :
		super( D, self).__init__( )
		self.volume = tuple( volume )
	def __or__( i, j ) :
		if not i.volume : return j
		return D( [ x | y
			for (x,y) in zip( i.volume, j.volume ) ] )
	def __le__( i, j ) :
		if not i.volume : return True
		if not j.volume : return False
		return all( [ x <= y
			for (x,y) in zip( i.volume, j.volume ) ] )
	def __eq__( i, j ) :
		return all( [ x == y
			for (x,y) in zip( i.volume, j.volume ) ] )
	def __ne__( i, j ) :
		return not ( i == j)
	def __repr__( self ) :
		if not self.volume : return "ø"
		return "x".join( [ str( x ) for x in self.volume ] )

## Implement an object-description mapping \delta: G\to (D,\sqcap)
## Basically an extesion of the dictonary to sets of keys
class Obj_Desc(object):
	def __init__( self, map ) :
		super(Obj_Desc, self).__init__()
		self.map = dict( map )
	def __call__( self, point = [] ) :
		if not point : return [D( )]
		return [ self.map[ k ]
			for k in self.map.keys( )
				if k in point ]
	def domain( self ) :
		return sorted( self.map.keys( ) )

## A calss for pattern structure (G, (D\sqcap), \delta)
class Pattern(object):
	def __init__(self, extent, intent, delta, name ):
		super(Pattern, self).__init__()
		self.intent = intent
		self.extent = extent
		self.name = name
		self.delta = delta
	def merge( this, that ) :
		extent = this.extent + that.extent
		map = dict( )
		for g in range( len( extent ) ) :
			e = extent[ g ]
			if e in this.extent :
				map[ g ] = this.delta.map[ this.extent.index( e ) ]
			else :
				map[ g ] = that.delta.map[ that.extent.index( e ) ]
		return Pattern( extent, this.intent, Obj_Desc( map ), this.name )
	def G(self, S):
		return [ self.extent[ x ] for x in S ]
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
## For any A ≤ G: uA = uduA; since A ≤ duA, AMON implies uduA ≤ uA whereas
##  by GALOIS from duA ≤ duA follows uA ≤ uduA. Similar result holds for
##  any description x: dx ≤ dudx.
	def u( pat, obj ) :
## Galois connection: join all descriptions
## The "up" function maps a subset of G to the minimally specific
##  description, through the map delta relation R:
##    uA = Sup{ delta(g) | g\in A }
		sup = pat.delta( [ ] )[ 0 ]
		for d in pat.delta( obj ) : sup |= d
		return sup
	def d( self, d ):
## The "down" map returns a subset of objects which fit the given
##  description d: dx = { g | delta(g) ≤ d }
		return set( [ g for g in self.delta.domain( )
			if self.delta( [ g ] )[ 0 ] <= d ] )
	def CbO( self ) :
## This function initiates the CbO for the Galois connection (u,d)
##  In this CbO on patterns extent and intent both grow simultaneously.
		return self.__CbO( set( [] ), D( [] ), self.delta.domain( ) )
	def lattice( self, reference = None ) :
## Construct a lattice, falsifying hypotheses if necessary
		return self.__lattice( set( [] ), D( [] ), {}, reference )
################################################################################
################################################################################
## Run the CbO
	def __CbO(self, A, B, G, lvl = 0 ):
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
## This is a recursive step of the CbO algorithm, it is called only when
##  a new concept has been canonically generated. Expect (A,B) to be a dual
##  pair with respect to the underlying Galois connection: B = uA, A = dB,
##  where u:Ω_A -> Ω_B is the Galois map to obtain the reflection of A in
##  the universe of B and d:Ω_B -> Ω_A reflects B into A's universe. G is
##  the size of A's universe (|Ω_A|)
##DEBUG##		assert( B == u(A) ) ; assert( A == d(B) )
##  Add the closure to the queue.
		Q = list( [ ( A, B ) ] )
## Examine elements in lexicographic order unless there is nothing more
##  to do. It is assumed that G, the set of objects, is sorted.
		while G :
			f = G.pop( 0 )
## There is no point in trying to expand by a member element, since we will generate
##  ourselves again
			if f in A : continue
## Close-by-One: generate a concept by adding a new member to it:
##   it is true that u(A | {f}) = uA & u{f} = B & u{f} (UNION)
			Z = B | self.u( [ f ] ) ; Y = self.d( Z )
			# Z = self.u( A | { f } ) ; Y = self.d( Z )
			# print " "*lvl, "%d»» %s %d]%s)" % ( lvl, "".join([str(e) for e in A]), f, "".join([str(g) for g in (Y-(A|{f}))]) )
## Perform the canonicity test: if all elements of Y \ A ( the newly generated
##  volume) are greater of equal to f then the generation (Y,Z) from A is
##  canonical. The test min( h | h in du(A | {f})-A ) ≥ f is actually equivalent
##  to checking if (du(A | {f})-A) & F is ø, where F={0..f-1}. In turn, since
##  C ≤ duC (CLOS) this test is equivalent to checking if A & F = du(A | {f}) & F
##  This insightful ``optimiation'' is possible only for objects in the pattern structure
			F = set( range( f ) )
			if F & Y == F & A :
## Every canonically generated pair (Y,uY) is guaranteed to have never been
##  generated by anyone before, because pairs are generated in lexicographic order.
##  Generate from the newly detected canonical dual pair (Y,uY), and add the proceeds
##  to the queue.
				Q += self.__CbO( Y, Z, list( G ), lvl + 1 )
		return Q
################################################################################
################################################################################
################################################################################
	def __lattice(self, A, B, L, cxt = None) :
################################################################################
## A simplified and tested version of [http://ami.hse.ru/data/165/315/1234/book.pdf] p. 117
################################################################################
## Create a new lattice node
		L[ label( A ) ] = lattice_node( label( A ), ( A, B ) )
## Collect all dual pairs that are generated from the current one by adding
##  a single new element 
		queue = dict( )
		for j in self.delta.domain( ) :
## Again, as in CbO above, there is no use in adding a part of the existing
##  dual pair
			if j in A : continue
## Like in CbO generate the next pair by adding a new element and closing with
##  the Galois connection (u,d). Note the join operator: it is here due to the nature
##  the description semi-lattice
			Y = B | self.u( [ j ] ) ; X = self.d( Y )
			# Y = self.u( A | { j } ) ; X = self.d( Y )
## Check if this closure is not falsified by the reference context
			if cxt and cxt.d( Y ) : continue
## Add a new heir to the list
			queue[ label( X ) ] = ( X, Y )
## Weed out nested heirs unfortunately due to the fact that set inclusion is 
##  inly a partial order on subsets, this has cubic complexity. This constructs
##  the next layer of nodes, which are related to the current via a covering
##  relation based on lattice partial order
		next = dict()
		for M in queue.values( ) :
## Scan the siblings for the smallest one with respect to the current M
			for Z in queue.values( ) :
				if Z[ 0 ] < M[ 0 ] : M = Z
## Add it to the call list.
			next[ label( M[0] ) ] = M
## Recursive call:
		for X,Y in next.values( ) :
## If the successor has already been visited, there is no need to traverse
##  its sub-lattice
			if label( X ) not in L :
				self.__lattice( X, Y, L, cxt )
## Connect this successor to the current node
			L[ label( A ) ].link_sup( L[ label( X ) ] )
		return L
	def save_lattice( self, L, filename ):
		if not L : return
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
			A, B = L[ np[ x ] ].pair
			lines.append( "%d '%s' '%s'\n" % ( x, ",".join( self.G( A ) ), B ) )
## Recover adjacency lists
		for k, x in L.items() :
			if x.sup :
				lines.append( "%d %s\n" %( pn[ k ], " ".join( [ str( pn[ z.tag ] ) for z in x.sup ] ) ) )
		with file( filename, "w+" ) as f:
			f.writelines( lines )

def skip_empty( f ):
# Skip an empty line if there is one
	pos = f.tell()
	if f.readline().rstrip():
		f.seek(pos)

def pattern_cxt(filename) :
	with file( filename, "r" ) as f:
# File format marker
		txt = f.readline().strip();
		if txt != "PATT":
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
## Read the data types of the intent: valid type are
##  S -- a finite set, I -- numeric inteval
		types = f.readline().strip()
# Load and process the object-attribute relationship,
#   expect that it is stored extent-wise
		cxt0 = dict()
		for g in range( G ) :
			descr = []
			input = zip( types, f.readline( ).split( "\t" ) )
			for d, v in input :
				if d is "S" :
					descr.append( S( v.split(",") ) )
				elif d is "I" :
					descr.append( I( [ float( v ) ] ) )
				else :
					raise Exception( "Unknown dscription type %s" % d )
## Instantiate the description
			cxt0[ g ] = D( descr )
		return Pattern( extent, intent, Obj_Desc( cxt0 ), name )

def usage() :
	print "%s <path-to-positive-context> [<path-to-positive-context>]" % (__file__)
	return -1

if __name__ == '__main__' :
	if len( sys.argv ) < 2 :
		exit( usage( ) )

	if len( sys.argv ) == 2 :
		filename = sys.argv[ 1 ]
		cxt = pattern_cxt( filename )
		cxt.save_lattice( cxt.lattice( ), filename + "-lattice.txt" )
		exit( 0 )

	filename_pos = sys.argv[ 1 ]
	filename_neg = sys.argv[ 2 ]
# filename_pos = "./ADACpos.cxt"
# filename_neg = "./ADACneg.cxt"
## Load the contexts
	pos = pattern_cxt( filename_pos )
	neg = pattern_cxt( filename_neg )

## Compute pattern concepts: positive context
	lat_pos = pos.lattice( )
	lat_neg = neg.lattice( )

## Closures are implications of the form d -> -/+
	pos.save_lattice( lat_pos, filename_pos + "-lattice.txt" )
	neg.save_lattice( lat_neg, filename_neg + "-lattice.txt" )

## Compute minimal unfalsified hypotheses
	hyp_pos = pos.lattice( neg )
	hyp_neg = neg.lattice( pos )

## Closures are implications of the form d -> -/+
	pos.save_lattice( hyp_pos, filename_pos + "-hypothses.txt" )
	neg.save_lattice( hyp_neg, filename_neg + "-hypothses.txt" )

## Join the contexts
	full = pos.merge( neg )
	full.save_lattice( full.lattice( ), filename_pos + "-full.txt" )

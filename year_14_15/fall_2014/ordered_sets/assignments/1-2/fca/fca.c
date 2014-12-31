#include <stdlib.h>
#include <stdio.h>
#include <string.h>


// g++ -O3 fca.c -o fca


// Define the context
#define M ( 9 )
#define G ( 8 )

const char CTX[G + 1][M + 1] = {
/*		    abcdefghi */
/* 1 */    "XX....X.."
/* 2 */  , "XX....XX."
/* 3 */  , "XXX...XX."
/* 4 */  , "X.X...XXX"
/* 5 */  , "XX.X.X..."
/* 6 */  , "XXXX.X..."
/* 7 */  , "X.XXX...."
/* 8 */  , "X.XX.X..."
/* 9 *//*, "XXXXXXXXX" */
};

/* -------------------------------------------------------------------------- *
 * --------------------------- Display primitives --------------------------- *
 * -------------------------------------------------------------------------- */
char *display_labels( char *buffer, const char *labels, const size_t n, const unsigned set )
{
	size_t j = 0;
	memset( buffer, '-', n );
	for( size_t i = 0; i < n; ++i )
		if( ( set >> i ) & 1 )
			buffer[ j++ ] = labels[ i ];
	if( !j ) ++j;
	buffer[ j ] = 0;
	return buffer;
}

inline void display_properties( const unsigned set, char *buffer )
{
	display_labels( buffer, "abcdefghijklmno", M, set );
}

inline void display_objects( const unsigned set, char *buffer )
{
	display_labels( buffer, "123456789ABCDEF", G, set );
}

/* -------------------------------------------------------------------------- *
 * ---------------------------- Poset primitives ---------------------------- *
 * -------------------------------------------------------------------------- */
inline bool subseteq( const unsigned a, const unsigned b )
{
// ({0,1}^n, subseteq) is a partialy ordered set with respect to the order
//  induced by set inclusion.
	return a == ( b & a );
}

/* -------------------------------------------------------------------------- *
 * ------- Galois connection between Posets of objects and attributes ------- *
 * -------------------------------------------------------------------------- */
// f^*:2^G -> 2^M is defined for any A subset of G as
//   f(A) = { m in M | object "g" has attribute "m" for all objects g in A }
unsigned properties( const unsigned objects )
{
	unsigned set = 0;
// Collect properties that are not possessed by at least one object
//  from the given set.
		for( int i = 0 ; i < G ; ++i ) {
// If the i-th object is not in the set, skip it.
			if( !( ( objects >> i ) & 1 ) ) continue;
// Collect all properties this object does not have
			unsigned mask = 0;
// To speed this step up just use a set of precalculated masks
			for( int j = 0 ; j < M ; ++j )
				mask |= ( CTX[i][j] == '.' ? 1 : 0 ) << j;
// ... accumulate properties not possessed by at least one object 
			set |= mask;
		}

// Flip the bits, since the the collected mask contains properties, which
//  are not possessed by at least one object.
	return ~set & ( ( 1 << M ) - 1 );
}

// f_*:2^M -> 2^G is defined for any B subset of M as
//   f(B) = { g in G | object "g" has attribute "m" for all attributes m in B }
unsigned objects( const unsigned properties )
{
	unsigned set = 0;
		for( int i = 0 ; i < M ; ++i ) {
			if( !( ( properties >> i ) & 1 ) ) continue;
			unsigned mask = 0;
// To speed this step up just use a set of precalculated TRNASPOSED masks
			for( int j = 0 ; j < G ; ++j )
				mask |= ( CTX[j][i] == '.' ? 1 : 0 ) << j;
			set |= mask;
		}

	return ~set & ( ( 1 << G ) - 1 );
}

// A closure operator on objects is the composition f_*(f^*(.))
inline unsigned object_closure( const unsigned set )
{
	return objects( properties( set ) );
}

// ... whereas a closure operator on attributes is the composition f^*(f_*(.))
inline unsigned property_closure( const unsigned set )
{
	return properties( objects( set ) );
}

// Lattice primitives
inline unsigned join( const unsigned a, const unsigned b )
{
	return b | a;
}

inline unsigned meet( const unsigned a, const unsigned b )
{
	return b & a;
}

/* -------------------------------------------------------------------------- *
 * --------------- Compute the formal concepts by brute force --------------- *
 * -------------------------------------------------------------------------- */
int main(int argc, char const *argv[])
{
	for(int i = argc ; i-- ; )
		printf( "%s\n", argv[ i ] );

	char buffer[2][16];
	FILE *fout;
	fout = fopen( "concepts.txt", "w+" );
// Build concepts. It does not matter whether we traverse the space of
//  objects or the space of properties.
	for( int i = 0 ; i < ( 1 << G ) ; ++i ) {
		unsigned p = properties( i );
		unsigned o = objects( p );
// The concept is a pair (g, m) such that g. = m and .m = g. Equivalently
//  a formal concept is a pair (g, g.) such that .(g.) = g, i.e g is closed
//  in the concept matrix. I wonder if one can define a topology on this
//  discrete poset.
		if( o != i ) continue;

		display_objects( i, buffer[ 0 ] );
		display_properties( p, buffer[ 1 ] );
		// display_objects( ( ~i ) & ( ( 1 << G ) - 1 ), buffer[ 0 ] );
		fprintf( fout, "%.*s\t%.*s\n", G, buffer[ 0 ], M, buffer[ 1 ] );
	}
	fclose( fout );

	return 0;
}

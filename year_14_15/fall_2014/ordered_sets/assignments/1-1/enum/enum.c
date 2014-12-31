#include <stdlib.h>
#include <stdio.h>
// #include <string.h>

/* FLAT bit array */
// Flat 1-D array primitives
//  assumes the size of the unsigned int is 32 bit
inline char get_bit( void * const data, const unsigned int i )
{
// Get the bin at position $\floor{ \frac{i}{ \text{sizeof(long)} * 8 / 1 } }$
// Fetch the bit at $i mod 32$;
	return ( ((unsigned int *)data)[ i >> 5 ] >> ( i & 31 ) ) & 1;
}

inline void set_bit( void *data, const unsigned int i )
{
	((unsigned int *)data)[ i >> 5 ] |= 1 << ( i & 31 );
}

inline void reset_bit( void *data, const unsigned int i )
{
	((unsigned int *)data)[ i >> 5 ] &= ~ ( 1 << ( i & 31 ) );
}

inline void toggle_bit( void *data, const unsigned int i )
{
	((unsigned int *)data)[ i >> 5 ] ^= ( 1 << ( i & 31 ) );
}

// Show 1-D array shaped like a NxN matrix
void show_bit( void *self, const unsigned int n )
{
// Print the header
	// printf( "#" );
	// for( int j = 0 ; j < 3 ; ++j ) printf( " %1d", j );
	// printf( "\n" );

	printf( " | " );
	for( int i = 0 ; i < n ; ++i ) {
		// printf( "%.1d", i );
		for( int j = 0 ; j < n ; ++j ) printf( "%c", get_bit( self, i + n * j )>0? '*' : '.' );
		printf( " | " );
	}
	printf( "\n" );
}

// Brute force check for reflexivity:
//  $\exists x\in A \text{s.t}  (x,x)\not\in R$
bool is_reflexive( void *const relation, const unsigned int n )
{
	for( int i = 0 ; i < n ; ++i )
		if( get_bit( relation, i * ( n + 1 ) ) == 0 ) return false;
	return true;
}

// Brute force check for antireflexivity:
//  $\exists x\in A s.t. (x,x)\in R$
bool is_antireflexive( void *const relation, const unsigned int n )
{
	for( int i = 0 ; i < n ; ++i )
		if( get_bit( relation, i * ( n + 1 ) ) != 0 ) return false;
	return true;
}

// Brute force check for symmetry:
//  $\exists x,y\in A s.t. (x,y)\not\in R or (y,x)\not\in R, or both$
bool is_symmetric( void *const relation, const unsigned int n )
{
	for( int i = 0 ; i < n ; ++i )
		for( int j = i+1 ; j < n ; ++j )
			if(    get_bit( relation, i * n + j ) != 0 
				|| get_bit( relation, j * n + i ) != 0 ) return false;
	return true;
}

// Brute force check for anit-symmetry:
//  $\exists x,y\in A, x\neq y s.t. (x,y)\in R and (y,x)\in R$
bool is_antisymmetric( void *const relation, const unsigned int n )
{
	for( int i = 0 ; i < n ; ++i )
		for( int j = i+1 ; j < n ; ++j )
			if(    get_bit( relation, i * n + j ) != 0 
				&& get_bit( relation, j * n + i ) != 0 ) return false;
	return true;
}

// Reduction of the check of asymmetry to two simlutaneous checks for
//  anti-reflexivity and anit-symmetry
bool is_asymmetric( void *const relation, const unsigned int n )
{
	return is_antireflexive( relation, n ) && is_antisymmetric( relation, n );
}

// Brute force check for transitivity:
// $\exists x,y,z\in A, s.t. (x,y), (y,z)\in R but (x,z)\not\in R$
bool is_transitive( void *const relation, const unsigned int n )
{
	char res = 0;
	for( int x = 0 ; x < n ; ++x )
		for( int y = 0 ; y < n ; ++y ) {
			if( get_bit( relation, x + y * n ) == 0 ) continue;
			for( int z = 0 ; z < n ; ++z ) {
				if( get_bit( relation, y + z * n ) == 0 ) continue;
				if( get_bit( relation, x + z * n ) == 0 ) return false;
			}
		}
	return true;
}

int main(int argc, char const *argv[])
{
// g++ -O3 ./enum.c -o enum
	for(int i = 0 ; i < argc ; ++i )
		printf( "%s\n", argv[ i ] );

// an unsigned long can hold up to an 8x8 binary matrix.
	unsigned long cnt = 0;

/* 25 bits would be enough to solve this task */
	for( unsigned long long a = 0 ; a < (1LL<<( 5 * 5 )) ; ++a ) {
		// if( !is_symmetric( &a, 5 ) || !is_reflexive( &a, 5 ) ) continue;
		if( !is_asymmetric( &a, 5 ) ) continue;
		if( !is_transitive( &a, 5 ) ) continue;
		show_bit( &a, 5 );
		++cnt;
	}
	
	printf("%ld\n", cnt);

	return 0;
} 


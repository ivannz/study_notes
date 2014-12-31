#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// Stringology #01-#02
// Knuth Morris Pratt algorithm: http://en.wikipedia.org/wiki/Knuth–Morris–Pratt_algorithm

struct _match_t
{
	size_t nMatches;
	size_t match[0];
};

typedef struct _match_t match_t;

size_t *boundary( char const *pattern, const size_t P )
{
// The alleged complexity of this algorithm is O(P)
	size_t *b = ( size_t * ) malloc( P * sizeof( size_t ) );

	for( size_t i = 0 ; i < P ; ++i ) {
// Set the length of the maximal bounding prefix of the subword[0..i]
//  to zero.
		b[ i ] = 0;

// Initialize the current symbol
		size_t k = i;
		while( k > 0 ) {
// The value b[ i ] points to the first symbol after the maximal bounding prefix
//  of the word P[ 0:i ]. In fact any bounding prefix of the sub-word formed by
//  a bounding prefix of P[ 0:i ] must itself be a bounding prefix of P[ 0:i ].
//  Indeed, if F = P[ 0:(f - 1) ] = P[ (i - (f - 1) ):i ] and S is a bounding
//  prefix of F then it must be true that S = F[ 0:s ] = F[ (f - 1 - s):(f - 1) ]
//  for some s. Since F[ 0:s ] = P[ 0:s ], F[ (f - 1 - s):(f - 1) ] must be
//  equal to P[ i - (f - 1) + (f - 1 - s):(f - 1) ] = P[ (i - s):i ]. Hence S
//  is a bounding prefix of P[0:i].
// Suppose F is the largest bounding prefix of P[0:i] and S is the largest
//  bounding prefix of F. If there A is another bounding prefix in P[0:i],
//  shorter than F, but longer than S, then it must be a bounding prefix of F,
//  whence S must be a bounding prefix of A, thus S is not maximal.

// Fetch the next greatest bounding prefix. If k = i then this gets the first
//  symbol after the maximal bounding prefix of the word P[0:(i-1)]. Otherwise
//  it gets the largest bounding prefix of the current bounding prefix.
			k = b[ k - 1 ];

// If the bounding prefix cannot be extended, then check the next one in order
			if( pattern[ i ] != pattern[ k ] ) continue;

//  ... otherwise the maximal bounding prefix of the word P[ 0:(k - 1) ] with
//  the new symbol P[ i ] is the maximal bounding prefix of the word P[ 0:i ].
			b[ i ] = k + 1;
			break;
		}
	}
	return b;
}

// Knuth Morris Pratt exact matching algorithm, which has complexity of O(T+P)
match_t *exact_match( char const *text, const size_t T, char const *pattern, const size_t P)
{
	match_t *matches = ( match_t * ) malloc( sizeof( match_t ) + T * sizeof( size_t ) );
	matches->nMatches = 0;

	if( T < P || P == 0 ) return matches;

// Compute the table of maximal bounding prefixes.
	size_t *b = boundary( pattern, P );

// Reset the current index in the text (i) and the current symbol (j) in the pattern.
	size_t i = 0, j = 0;

// If there is not enough symbols left in the text for a pattern to match,
//  then there is no point in continuing the search. Otherwise go on until
//  the currently matched symbol in the text is actually within the text.
	while( i + j < T && i + P < T + 1 ) {
// Scan until the first mismatch or the end of pattern
		while( j < P && ( text[ i + j ] == pattern[ j ] ) ) j++;

// The match is successful, only if the whole pattern has been scanned.
		if( j == P )
			matches->match[ matches->nMatches++ ] = i;

// Partial prefix matches are actually no different form the whole pattern
//  match. Just find the maximal bounding prefix of the matched part of the
//  pattern, and shift the position in the text forward, so that the boundary
//  prefix at the head of the pattern is aligned with the boundary in the
//  matched part of the text.
		printf("%d: %.*s (%c) -> ", j, P, text + i, text[ i+j ] );
		if( j > 0 ) {
			i += j - b[ j - 1 ];
// Shift the current symbol of the pattern back to the first symbol after the
//  maximal boundary of the part matched so far.
			j = b[ j - 1 ];
		} else {
// Failed to match anything. Move to the next symbol in the text and restart the
//  matching of the pattern.
			i++;
			j = 0;
		}
		printf("%d: %.*s (%c) @ %d\n", j, P, text + i, text[ i+j ], i+j );
	}

	free( b );

	return matches;
}

// ./lect01 abacabadabacabaeabacabadabacabafabacabadabacabaeabacabadabacabagabacabadabacabaeabacabadabacabafabacabadabacabaeabacabadabacaba llolkloll
// ./lect01 abacabadabacabaeabacabadabacabafabacabadabacabaeabacabadabacabagabacabadabacabaeabacabadabacabafabacabadabacabaeabacabadabacaba abadaba
// gcc -O3 lect01.c -o lect01
int main(int argc, char const *argv[])
{
	for( int i = argc ; i-- ; ) {
		printf( "%s\n", argv[ i ] );
		size_t *b = boundary( argv[ i ], strlen( argv[ i ] ) );

		for( int j = 0 ; j < strlen( argv[ i ] ) ; ++j )
			printf("%2d ", b[ j ] );
		printf( "\n" );

		free( b );
	}

		size_t *b = boundary( "participate in parachute", 24 );

		for( int j = 0 ; j < 24 ; ++j )
			printf("%c %2d\n", "participate in parachute"[ j ], b[ j ] );
		printf( "\n" );
		free( b );

// Do some matching
	match_t *matches = exact_match( argv[ 1 ], strlen( argv[ 1 ] ), argv[ 2 ], strlen( argv[ 2 ] ) );

	for( size_t m = 0 ; m < matches->nMatches ; ++m ) {
		printf("%d\n", matches->match[ m ] );
	}

	free( matches );
	return 0;
}
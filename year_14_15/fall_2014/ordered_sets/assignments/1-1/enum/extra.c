/* bit_matrix */
struct _bit_matrix_t {
	unsigned int rows, cols;
	char data[0];
};
typedef struct _bit_matrix_t bit_matrix;

bit_matrix * bit_matrix_alloc(const unsigned int rows, const unsigned int cols)
{
	bit_matrix *self = (bit_matrix *) malloc( sizeof( bit_matrix ) + rows * cols * sizeof( char ) );
	self->cols = cols ; self->rows = rows ;

	memset( self->data, -1, rows * cols * sizeof( char ) );

	return self;
}

void bit_matrix_dealloc(bit_matrix * self)
{
	if( self == NULL ) return;
	self->cols = self->rows = 0;
	free( self );
}

char bit_matrix_get( bit_matrix * const self, unsigned int i, unsigned int j )
{
// The matri elelemts are stored in the array column-wise. Thus in order to access
//  the lemenet at the i-th row and j-th column one should actualy access
//  j*rows + i
	if( i < self->rows && j < self->cols )
		return self -> data[ self -> rows * j + i ];
	return -1;
}

void bit_matrix_set( bit_matrix * self, unsigned int i, unsigned int j, char value )
{
// The matri elelemts are stored in the array column-wise. Thus in order to access
//  the lemenet at the i-th row and j-th column one should actualy access
//  j*rows + i
	if( i < self->rows && j < self->cols )
		self->data[ self -> rows * j + i ] = value;
}

void bit_matrix_show( bit_matrix * const self)
{
// Print the header
	printf( "##" );
	for( int j = 0 ; j < self->cols ; ++j ) printf( "\t%2d", j );
	printf( "\n" );

	for( int i = 0 ; i < self->rows ; ++i ) {
		printf( "%4d", i );
		for( int j = 0 ; j < self->cols ; ++j ) printf( "\t%-+4d", bit_matrix_get( self, i, j ) );
		printf( "\n" );
	}
}


bool is_reflexive_const( void *const relation, const unsigned int n )
{
	char res = 0;
	for( int i = 0 ; i < n ; ++i )
		res |= ( 1 ^ get_bit( relation, i * ( n + 1 ) ) );
	return res == 0;
}
bool is_antireflexive_const( void *const relation, const unsigned int n )
{
	char res = 0;
	for( int i = 0 ; i < n ; ++i )
		res |= ( 0 ^ get_bit( relation, i * ( n + 1 ) ) );
	return res == 0;
}

bool is_symmetric_const( void *const relation, const unsigned int n )
{
	char res = 0;
	for( int i = 0 ; i < n ; ++i )
		for( int j = i+1 ; j < n ; ++j )
			res |= get_bit( relation, i * n + j ) ^ get_bit( relation, j * n + i );

	return res == 0;
}

bool is_antisymmetric_const( void *const relation, const unsigned int n )
{
	char res = 0;
	for( int i = 0 ; i < n ; ++i )
		for( int j = i+1 ; j < n ; ++j )
// min r_{x,y}, r_{y,x} = 0
			res |= ( get_bit( relation, i * n + j ) & get_bit( relation, j * n + i ) );

	return res == 0;
}

bool is_asymmetric_const( void *const relation, const unsigned int n )
{
	return is_antireflexive_const( relation, n ) && is_antisymmetric_const( relation, n );
}

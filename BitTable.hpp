/**
  Use an 64bit integer to represent bit table.

  Li Song
  July 24, 2012

  - Modified again from May 27, 2017.
*/


#ifndef _BIT_TABLE_HEADER
#define _BIT_TABLE_HEADER

#include <stdio.h>

typedef unsigned long long int UINT64 ;
#define UNIT_SIZE (sizeof( UINT64 ) * 8 )
#define UNIT_MASK ((UINT64)63) 

class BitTable
{
private:
        UINT64 *tab ; // The bits
        int size ; // The size of the table
	int asize ; // The size of the array (size/64).
public:
        BitTable() // Intialize a bit table.
	{
		//tab = new UINT64[1] ;
		//size = UNIT_SIZE ;
		tab = NULL ;
		size = 0 ;
	}

        BitTable( int s ) // Initalize a bit table with size s.
	{
		if ( s == 0 )
			s = UNIT_SIZE ;

		if ( s & UNIT_MASK )
		{
			asize = s / UNIT_SIZE + 1 ;
			tab = new UINT64[ asize ]() ;
		}
		else
		{
			asize = s / UNIT_SIZE ;
			tab = new UINT64[ asize ]() ;
		}

		size = s ;
		Reset() ;
	}

        ~BitTable() 
	{
		//printf( "hi %d\n", size ) ;
		//printf( "%s\n", __func__ ) ;
		//if ( tab != NULL )
		//	delete [] tab ;
		//size = -1 ;
		
		// Manually release the memory.
	}
	
	void Init( int s ) // Initialize a bit table with size s.
	{
		if ( tab != NULL )
			delete [] tab ;

		if ( s == 0 )
			s = UNIT_SIZE ;

		if ( s & UNIT_MASK )
		{
			asize = s / UNIT_SIZE + 1 ;
			tab = new UINT64[ asize ]() ;
		}
		else
		{
			asize = s / UNIT_SIZE ;
			tab = new UINT64[ asize ]() ;
		}

		size = s ;
		Reset() ;
	}

	void Release()
	{
		if ( tab != NULL )
			delete[] tab ;
	}

	void Reset()  // Make every value 0.
	{
		for ( int i = 0 ; i < asize ; ++i )
			tab[i] = 0 ;
	}

        void Set( int i ) // Set the ith bit.
	{
		int ind, offset ;
		ind = i / UNIT_SIZE ;
		offset = i & UNIT_MASK ;
		//printf( "%d,%d %d,%d: %llx", 311 / UNIT_SIZE, 311 & UNIT_MASK, 279 / UNIT_SIZE, 279 & UNIT_MASK, tab[ind] ) ; 	
		tab[ind] |= ( (UINT64)1 << offset ) ;
		//printf( " %llx %d\n", tab[ind], UNIT_SIZE ) ;
	}

        void Unset( int i ) // Unset the ith bit
	{
		int ind, offset ;
		ind = i / UNIT_SIZE ;
		offset = i & UNIT_MASK ;

		tab[ind] &= ( (UINT64)-1 ^ ( (UINT64)1 << offset ) ) ;
	}

        void Flip( int i ) // Flip the ith bit. Same as xor 1.        
	{
		int ind, offset ;
		ind = i / UNIT_SIZE ;
		offset = i & UNIT_MASK ; 	

		tab[ind] ^= ( (UINT64)1 << offset ) ;
	}

        bool Test( int i ) const // Test the ith bit.
	{
		if ( i >= size )
			return false ;

		int ind, offset ;
		ind = i / UNIT_SIZE ;
		offset = i & UNIT_MASK ; 

		if ( tab[ind] & ( (UINT64)1 << offset ) )
			return true ;
		else
			return false ;
	}

	void Not() // Not all the bits.
	{
		int i ;
		for ( i = 0 ; i < asize ; ++i )
		{
			tab[i] = ~tab[i] ;
		}
	}

	void And( const BitTable &in ) // Do the "and" on each bits.
	{
		if ( asize != in.asize )
			return ;
		int i ;
		for ( i = 0 ; i < asize ; ++i )
			tab[i] &= in.tab[i] ; 
	}

	void Or( const BitTable &in ) // Do the "or" on each bits.
	{
		if ( asize != in.asize )
			return ;
		int i ;
		for ( i = 0 ; i < asize ; ++i )
			tab[i] |= in.tab[i] ; 
	}

	void Xor( const BitTable &in )
	{
		if ( asize != in.asize )
			return ;
		int i ;
		for ( i = 0 ; i < asize ; ++i )
			tab[i] ^= in.tab[i] ;
	}
	
	// Unset all the bits outside [s,e]
	void MaskRegionOutside( int s, int e )
	{
		int i ;
		// mask out [0, s-1].
		int ind, offset ;
		if ( s > 0 )
		{
			ind = (s - 1) / UNIT_SIZE ;
			offset = ( s - 1 ) & UNIT_MASK ;
			for ( i = 0 ; i < ind - 1 ; ++i )			
				tab[i] = 0 ;
			tab[i] = ( tab[i] >> ( offset + 1 ) ) << ( offset + 1 ) ; 
		}
		
		// mask out [e+1, size-1]
		if ( e < size - 1 )
		{
			ind = ( e + 1 ) / UNIT_SIZE ;
			offset = ( e + 1 ) & UNIT_MASK ;
			for ( i = ind + 1 ; i < asize ; ++i )
				tab[i] = 0 ;
			tab[ind] = ( tab[ind] << ( UNIT_SIZE - offset ) ) >> ( UNIT_SIZE - offset ) ;
		}
	}

	void ShiftOneLeft()
	{
		int i ;
		UINT64 carry = 0 ;
		UINT64 lastTabMask ;
		if ( ( size & UNIT_MASK ) == 0 )
			lastTabMask = (UINT64)(-1) ;
		else
			lastTabMask = ( 1 << ( size & UNIT_MASK ) ) - 1 ; 

		for ( i = 0 ; i < asize - 1 ; ++i )
		{
			UINT64 tmp = ( tab[i] >> ( UNIT_SIZE - 1 ) ) & 1 ;
			tab[i] = ( tab[i] << 1 ) | carry ;
			carry = tmp ;
		}
		tab[i] = ( ( tab[i] << 1 ) | carry ) & lastTabMask ;
	}

	void ShiftOneRight()
	{
		int i ;
		UINT64 carry = ( tab[ asize - 1 ] & 1 ) ;
		tab[ asize - 1 ] >>= 1 ;
		for ( i = asize - 2 ; i >= 0 ; --i )
		{
			UINT64 tmp = ( tab[i] & 1 ) ;
			tab[i] = ( tab[i] >> 1 ) | ( carry << ( UNIT_SIZE - 1 ) ) ;
			carry = tmp ;
		}
	}

	bool IsAllZero()
	{
		int i ;
		for ( i = 0 ; i < asize ; ++i )
			if ( tab[i] != 0 )
				return false ;
		return true ;
	}

	int Count() const // Count how many 1.
	{
		if ( size <= 0 )
			return 0 ;
		UINT64 k ;
		int i, ret = 0 ;
		for ( i = 0 ; i < asize - 1 ; ++i )
		{
			k = tab[i] ;
			while ( k )
				//for ( j = 0 ; j < UNIT_SIZE ; ++j ) 	
			{
				if ( k & 1 )
					++ret ;
				//printf( "### %d %d %d %d\n", ret, asize, size, k ) ;
				k /= 2 ;
			}
		}
		//printf( "(%d) ", ret ) ;	
		k = tab[ asize - 1 ] ;
		for ( i = 0 ; i < (int)( size & UNIT_MASK ) ; ++i )
		{
			if ( k & 1 )
				++ret ;
			k /= 2 ;
		}
		//for ( i = 0 ; i < asize ; ++i ) 
		//	printf( "%llu ", tab[i] ) ;
		//printf( "(%d)\n", ret ) ;
		return ret ;
	}

	void GetOnesIndices( std::vector<int> &indices )
	{
		if ( size <= 0 )
			return ;
		UINT64 k ;
		int i ;
		int ind = 0 ;
		for ( i = 0 ; i < asize - 1 ; ++i )
		{
			k = tab[i] ;
			ind = i * UNIT_SIZE ; 
			while ( k )
				//for ( j = 0 ; j < UNIT_SIZE ; ++j ) 	
			{
				if ( k & 1 )
					indices.push_back( ind ) ;
				//printf( "### %d %d %d %d\n", ret, asize, size, k ) ;
				k /= 2 ;
				++ind ;
			}
		}
		//printf( "(%d) ", ret ) ;	
		k = tab[ asize - 1 ] ;
		ind = i * UNIT_SIZE ; 
		for ( i = 0 ; i < (int)( size & UNIT_MASK ) ; ++i )
		{
			if ( k & 1 )
				indices.push_back( ind ) ;
			k /= 2 ;
			++ind ;
		}
	}

        bool IsEqual( const BitTable &in ) const// Test wether two bit tables equal.   
	{
		int i ;
		if ( in.size != size )
			return false ;
		//printf( "== %d %d\n", in.size, size ) ;
		//if ( size & UNIT_MASK )
		//	k = size / UNIT_SIZE + 1 ;
		//else
		//	k = size / UNIT_SIZE ;	

		for ( i = 0 ; i < asize ; ++i )	
			if ( tab[i] != in.tab[i] )
				return false ;
		return true ;
	}
	
	// Return the location of the first difference. -1 if the same.
	int GetFirstDifference( const BitTable &in ) const
	{
		int as = asize < in.asize ? asize : in.asize ;
		//int s = size < in.size ? size : in.size ;
		int i, j ;

		for  ( i = 0 ; i < as ; ++i )
		{
			if ( tab[i] != in.tab[i] )
			{
				int k1 = tab[i] ;					
				int k2 = in.tab[i] ;
				for ( j = 0 ; j < (int)UNIT_SIZE ; ++j )
				{
					if ( ( k1 & 1 ) != ( k2 & 1 ) )
					{
						return j + UNIT_SIZE * i ;
					}
					k1 >>= 1 ;
					k2 >>= 1 ;
				}
			}
		}
		
		for ( ; i < asize ; ++i )
		{
			if ( tab[i] != 0 )
			{
				int k = tab[i] ;
				for ( j = 0 ; j < (int)UNIT_SIZE ; ++j )
				{
					if ( k & 1 )
						return j + UNIT_SIZE * i ;
					k >>= 1 ;
				}
			}
		}
		
		for ( ; i < in.asize ; ++i )
		{
			if ( in.tab[i] != 0 )
			{
				int k = in.tab[i] ;
				for ( j = 0 ; j < (int)UNIT_SIZE ; ++j )
				{
					if ( k & 1 )
						return j + UNIT_SIZE * i ;
					k >>= 1 ;
				}
			}
		}
		return -1 ;
	}

	void Duplicate( BitTable &in )
	{
		if ( tab != NULL )
			delete [] tab ;
		int i ;
		size = in.size ;
		asize = in.asize ;
		tab = new UINT64[ asize ]() ;
		for ( i = 0 ; i < asize ; ++i )
			tab[i] = in.tab[i] ;
	}

	void Assign( BitTable &in )
	{
		int i ;
		for ( i = 0 ; i < asize ; ++i )
			tab[i] = in.tab[i] ;
	}
} ;


#endif

#include "Constraints.hpp"

// return whether this constraint is compatible with the subexons.
bool Constraints::ConvertAlignmentToBitTable( struct _pair *segments, int segCnt, 
	struct _subexon *subexons, int seCnt, int seStart, struct _constraint &ct ) 
{
	int i, k ;
	k = 0 ;
	ct.vector.Init( seCnt ) ;
	// Each segment of an alignment can cover several subexons.
	// But the first and last segment can partially cover a subexon.
	for ( i = 0 ; i < segCnt ; ++i )
	{
		int leftIdx, rightIdx ; // the range of subexons covered by this segment.
		leftIdx = -1 ;
		rightIdx = -1 ;
		
		for ( ; k < seCnt ; ++k )
		{
			//printf( "(%d:%d %d):(%d:%d %d)\n", i, segments[i].a, segments[i].b, k, subexons[k].start, subexons[k].end ) ;
			if ( subexons[k].start > segments[i].b )
				break ;
			
			if ( ( subexons[k].start >= segments[i].a && subexons[k].end <= segments[i].b ) 
				|| ( i == 0 && subexons[k].start < segments[i].a && subexons[k].end <= segments[i].b ) 
				|| ( i == segCnt - 1 && subexons[k].start >= segments[i].a && subexons[k].end > segments[i].b ) )
			{
				if ( leftIdx == -1 )
					leftIdx = k ;
				rightIdx = k ;
				ct.vector.Set( k ) ;
			}
			else
			{
				return false ;
			}
		}

		if ( leftIdx == -1 )
			return false ;

		// The cover contradict the boundary.
		if ( !( ( subexons[leftIdx].leftType == 0 || subexons[leftIdx].start <= segments[i].a )
			&& ( subexons[rightIdx].rightType == 0 || subexons[rightIdx].end >= segments[i].b ) ) )
			return false ;
	}

	return true ;
}

void Constraints::CoalesceSameConstraints()
{
	int i, k ;
	std::sort( constraints.begin(), constraints.end(), CompSortConstraints ) ;
	int size = constraints.size() ;
	k = 0 ;
	for ( i = 1 ; i < size ; ++i )
	{
		if ( constraints[k].vector.IsEqual( constraints[i].vector ) )
		{
			constraints[k].support += constraints[i].support ; 
			constraints[i].vector.Release() ;
		}
		else
		{
			++k ;
			if ( k != i )
				constraints[k] = constraints[i] ;
		}
	}
	constraints.resize( k + 1 ) ;
}

int Constraints::BuildConstraints( struct _subexon *subexons, int seCnt, int start, int end )
{
	int i ;
	int tag = 0 ;

	// 
	int size = constraints.size() ;
	if ( size > 0 )
	{
		for ( i = 0 ; i < size ; ++i )
			constraints[i].vector.Release() ;
		std::vector<struct _constraint>().swap( constraints ) ;
	}

	while ( alignments.Next() != -1 )
	{
		if ( alignments.GetChromId() < subexons[0].chrId )
			continue ;
		else if ( alignments.GetChromId() > subexons[0].chrId )
			break ;

		// locate the first subexon in this region that overlapps with current alignment.
		for ( ; tag < seCnt && subexons[tag].end < alignments.segments[0].a ; ++tag )
			;
		
		if ( tag >= seCnt )
			break ;
		
		if ( alignments.segments[ alignments.segCnt - 1 ].b < subexons[tag].start )
			continue ;

		struct _constraint ct ;
		ct.vector.Init( seCnt ) ;
		//printf( "%s %d: %d-%d | %d-%d\n", __func__, alignments.segCnt, alignments.segments[0].a, alignments.segments[0].b, subexons[tag].start, subexons[tag].end ) ;
		if ( ConvertAlignmentToBitTable( alignments.segments, alignments.segCnt, 
				subexons, seCnt, tag, ct ) )
		{	
			constraints.push_back( ct ) ; // if we just coalesced but the list size does not decrease, this will force capacity increase.
			// TODO: what happens it does not coalesce enough?
			if ( constraints.size() == constraints.capacity() )
			{
				//printf( "start coalescing.\n" ) ;
				CoalesceSameConstraints() ;	
			}
		}
		else
			ct.vector.Release() ;
		
		// TODO: add the mate-pair information.
	}
	return 0 ;
}

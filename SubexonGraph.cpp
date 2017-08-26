#include "SubexonGraph.hpp"

void SubexonGraph::GetGeneBoundary( int tag, int strand, int &boundary, int timeStamp )
{
	if ( visit[tag] == timeStamp )	
		return ;
	//printf( "%d %d\n", tag, timeStamp ) ;
	visit[tag] = timeStamp ;
	if ( subexons[tag].end > boundary )
		boundary = subexons[tag].end ;
	//if ( subexons[tag].start == 2858011 )
	//	printf( "%d: %d %d\n", tag, subexons[tag].nextCnt, subexons[tag].prevCnt) ;
	int i ;
	if ( IsSameStrand( subexons[tag].rightStrand, strand ) )
	{
		int cnt = subexons[tag].nextCnt ;
		for ( i = 0 ; i < cnt ; ++i )
		{
			//printf( "next of %d: %d %d\n", tag, i, subexons[tag].next[i] ) ;
			GetGeneBoundary( subexons[tag].next[i], strand, boundary, timeStamp ) ;
		}
	}
	else
	{
		// The only way to reach this subexon is through the adjacent exon bouldary
		// like ]...]. So we can regard this subexon as unvisited, and let the
		// procecure in GetGenintervalIdx to decide whether we want to include
		// this subexon or not.
		visit[tag] = -1 ;
	}

	if ( IsSameStrand( subexons[tag].leftStrand, strand ) )
	{
		int cnt = subexons[tag].prevCnt ;
		for ( i = 0 ; i < cnt ; ++i )
			GetGeneBoundary( subexons[tag].prev[i], strand, boundary, timeStamp ) ;
	}
	else
	{
		visit[tag] = -1 ;
	}
}

int SubexonGraph::GetGeneIntervalIdx( int startIdx, int &endIdx, int timeStamp )
{
	int i ;
	int seCnt = subexons.size() ;
	if ( startIdx >= seCnt )
		return -1 ;
	int farthest = -1 ;
	GetGeneBoundary( startIdx, subexons[ startIdx ].rightStrand, farthest, timeStamp ) ;

	for ( i = startIdx ; i < seCnt ; ++i )
	{
		if ( subexons[i].start > farthest || subexons[i].chrId != subexons[ startIdx ].chrId )								
			break ;
		GetGeneBoundary( i, subexons[i].rightStrand, farthest, timeStamp ) ;
	}
	endIdx = i - 1 ;

	return endIdx ;
}

int SubexonGraph::ComputeGeneIntervals()
{
	int i, cnt ;
	int seCnt = subexons.size() ;
	visit = new int[seCnt] ;
	memset( visit, -1, sizeof( int ) * seCnt ) ;
	int tag = 0 ;
	cnt = 0 ;
	while ( 1 )		
	{
		struct _geneInterval ngi ;
		if ( GetGeneIntervalIdx( tag, ngi.endIdx, cnt ) == -1 )
			break ;
		++cnt ;
		ngi.startIdx = tag ;
		ngi.start = subexons[ ngi.startIdx ].start ;
		ngi.end = subexons[ ngi.endIdx ].end ;
		
		tag = ngi.endIdx + 1 ;
		// Adjust the extent
		// Adjust the start 
		if ( subexons[ ngi.startIdx ].leftStrand != 0 
			&& subexons[ngi.startIdx].leftStrand != subexons[ngi.startIdx ].rightStrand )
			// We should make sure that rightstrand is non-zero whenever left-strand is non-zero for the startIdx.
		{
			for ( i = ngi.startIdx ; i >= 0 ; --i )	
			{
				if ( ( subexons[i].leftType == 1 && subexons[i].leftClassifier < classifierThreshold ) // an end within the subexon
					|| ( subexons[i].leftType == 0 ) // probably a overhang subexon. It should be a subset of the criterion following.
				  	|| ( i > 0 && subexons[i - 1].end + 1 < subexons[i].start ) ) // a gap. 
					break ;
			}
			ngi.start = subexons[i].start ;
		}

		// Adjust the end. And here, we also need to decide wether we need to adjust "tag" or not, 
		// because the next interval might be overlap with current interval by the last subexon.
		if ( subexons[ ngi.endIdx ].rightStrand != 0 
			&& subexons[ngi.endIdx].leftStrand != subexons[ngi.endIdx ].rightStrand )
		{
			for ( i = ngi.endIdx ; i < seCnt ; ++i )	
			{
				if ( ( subexons[i].rightType == 2 && subexons[i].rightClassifier < classifierThreshold ) // an end within the subexon
					|| ( subexons[i].rightType == 0 ) // probably a overhang subexon.
				  	|| ( i < seCnt - 1 && subexons[i].end + 1 < subexons[i + 1].start ) ) // a gap
					break ;
			}
			ngi.end = subexons[i].end ;
			
			if ( subexons[ ngi.endIdx ].rightType == 2 )
			{
				for ( i = ngi.endIdx ; i >= ngi.startIdx ; --i )
				{	
					if ( subexons[i].leftType == 1 )
						break ;
				}
				// The last region overlapps.
				if ( i >= ngi.startIdx && subexons[i].leftStrand != subexons[ ngi.endIdx ].rightStrand )
					--tag ;
			}
		}
		geneIntervals.push_back( ngi ) ;
	}
	delete[] visit ;

	return cnt ;
}

int SubexonGraph::ExtractSubexons( int startIdx, int endIdx, struct _subexon *retList )
{
	int i, j, k ;
	int cnt = endIdx - startIdx + 1 ;
	//printf( "%s: %d %d %d\n", __func__, startIdx, endIdx, cnt ) ;
	for ( i = 0 ; i < cnt ; ++i )
	{
		retList[i] = subexons[i + startIdx] ;
		retList[i].prev = new int[ retList[i].prevCnt ]	;
		retList[i].next = new int[ retList[i].nextCnt ] ;
		
		for ( j = 0 ; j < retList[i].prevCnt ; ++j )
			retList[i].prev[j] = subexons[i + startIdx].prev[j] - startIdx ;
		for ( j = 0 ; j < retList[i].nextCnt ; ++j )
			retList[i].next[j] = subexons[i + startIdx].next[j] - startIdx ;
		
		for ( j = 0, k = 0 ; j < retList[i].prevCnt ; ++j )
			if ( retList[i].prev[j] >= 0 && retList[i].prev[j] < cnt )
			{
				retList[i].prev[k] = retList[i].prev[j] ;
				++k ;
			}
		retList[i].prevCnt = k ;

		for ( j = 0, k = 0 ; j < retList[i].nextCnt ; ++j )
			if ( retList[i].next[j] >= 0 && retList[i].next[j] < cnt )
			{
				retList[i].next[k] = retList[i].next[j] ;
				++k ;
			}
		retList[i].nextCnt = k ;
	}
	return cnt ;	
}

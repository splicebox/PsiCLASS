#include "TranscriptDecider.hpp"

#define USE_DP 200000

void TranscriptDecider::OutputTranscript( FILE *fp, int baseGeneId, struct _subexon *subexons, struct _transcript &transcript )
{
	int i, j ;
	// determine the strand
	std::vector<int> subexonInd ;
	transcript.seVector.GetOnesIndices( subexonInd ) ;

	// Determine the strand
	char strand[2] = "." ;
	int size = subexonInd.size() ;
	if ( size > 1 )
	{
		// locate the intron showed up in this transcript.
		for ( i = 0 ; i < size - 1 ; ++i )
		{
			int nextCnt = subexons[ subexonInd[i] ].nextCnt ;
			if ( nextCnt == 0 )
				continue ;

			for ( j = 0 ; j < nextCnt ; ++j )
				if ( transcript.seVector.Test( subexons[ subexonInd[i] ].next[j] ) )
					break ;
			if ( j < nextCnt )
			{
				if ( subexons[ subexonInd[i] ].rightStrand == 1 )
					strand[0] = '+' ;
				else if (  subexons[ subexonInd[i] ].rightStrand == -1 )
					strand[0] = '-' ;
				break ;
			}
		}
	}

	// TODO: transcript_id
	char *chrom = alignments.GetChromName( subexons[0].chrId ) ;
	int a = subexonInd[0] ;
	int b = subexonInd[size - 1] ;
	char prefix[10] = "" ;

	fprintf( fp, "%s\tCLASSES\ttranscript\t%d\t%d\t1000\t%s\t.\tgene_id \"%s%s.%d\"; transcript_id \"%s%s.%d.%d\"; Abundance \"%.6lf\";\n",
			chrom, subexons[a].start + 1, subexons[b].end + 1, strand,
			prefix, chrom, geneId[a],
			prefix, chrom, geneId[a], transcriptId[ geneId[a] - baseGeneId ], transcript.abundance ) ;
	for ( i = 0 ; i < size ; ++i )
	{
		fprintf( fp, "%s\tCLASSES\texon\t%d\t%d\t1000\t%s\t.\tgene_id \"%s%s.%d\"; "
				"transcript_id \"%s%s.%d.%d\"; exon_number \"%d\"; Abundance \"%.6lf\"\n",
				chrom, subexons[ subexonInd[i] ].start + 1, subexons[ subexonInd[i] ].end + 1, strand,
				prefix, chrom, geneId[a],
				prefix, chrom, geneId[a], transcriptId[ geneId[a] - baseGeneId ],
				i + 1, transcript.abundance ) ;
	}
	++transcriptId[ geneId[a] - baseGeneId ] ;
}

void TranscriptDecider::SetGeneId( int tag, struct _subexon *subexons, int id )
{
	if ( geneId[tag] != -1 )
		return ;
	int i ;
	geneId[ tag ] = id ;
	int cnt = subexons[tag].nextCnt ;
	for ( i = 0 ; i < cnt ; ++i )
		SetGeneId( subexons[tag].next[i], subexons, id ) ;

	cnt = subexons[tag].prevCnt ;
	for ( i = 0 ; i < cnt ; ++i )
		SetGeneId( subexons[tag].prev[i], subexons, id ) ;
}

void TranscriptDecider::InitTranscriptId( int baseGeneId, int usedGeneId )
{
	int i ;
	for ( i = 0 ; i < usedGeneId - baseGeneId ; ++i )
		transcriptId[i] = 0 ;
}

// Return 0 - uncompatible or does not overlap at all. 1 - fully compatible. 2 - Head of the constraints compatible with the tail of the transcript
// the partial compatible case (return 2) mostly likely happen in DP where we have partial transcript.
int TranscriptDecider::IsConstraintInTranscript( struct _transcript transcript, struct _constraint &c ) 
{
	//printf( "%d %d, %d %d\n", c.first, c.last, transcript.first, transcript.last ) ;
	if ( c.first < transcript.first || c.first > transcript.last ) // no overlap or starts too early.
		return 0 ; 

	// Extract the subexons we should focus on.
	int s, e ;
	s = c.first ;
	e = c.last ;
	bool returnPartial = false ;
	if ( e > transcript.last ) // constraints ends after the transcript.
	{
		if ( transcript.partial )	
		{
			e = transcript.last ;
			returnPartial = true ;
		}
		else
			return 0 ;
	}
	/*printf( "%s: %d %d: (%d %d) (%d %d)\n", __func__, s, e,
	  transcript.seVector.Test(0), transcript.seVector.Test(1), 
	  c.vector.Test(0), c.vector.Test(1) ) ;*/

	compatibleTestVectorT.Assign( transcript.seVector ) ;
	compatibleTestVectorT.MaskRegionOutside( s, e ) ;

	compatibleTestVectorC.Assign( c.vector ) ;
	if ( e > transcript.last )
		compatibleTestVectorC.MaskRegionOutside( s, e ) ;
	/*printf( "after masking: (%d %d) (%d %d)\n", 
	  compatibleTestVectorT.Test(0), compatibleTestVectorT.Test(1), 
	  compatibleTestVectorC.Test(0), compatibleTestVectorC.Test(1) ) ;*/

	// Test compatible.
	int ret = 0 ;
	if ( compatibleTestVectorT.IsEqual( compatibleTestVectorC ) )
	{
		if ( returnPartial )
			ret = 2 ;
		else
			ret = 1 ;
	}

	return ret ;
}

int TranscriptDecider::SubTranscriptCount( int tag, struct _subexon *subexons, int *f )
{
	if ( f[tag] != -1 )
		return f[tag] ;

	int ret = 0 ;
	int i ;
	if ( subexons[tag].canBeEnd )
		ret = 1 ;
	for ( i = 0 ; i < subexons[tag].nextCnt ; ++i )
	{
		ret += SubTranscriptCount( subexons[tag].next[i], subexons, f ) ;
	}

	if ( ret == 0 )
		ret = 1 ;
	return f[tag] = ret ;
}

void TranscriptDecider::CoalesceSameTranscripts( std::vector<struct _transcript> &t ) 
{
	int i, k ;
	std::sort( t.begin(), t.end(), CompSortTranscripts ) ;

	int size = t.size() ;
	k = 0 ;
	for ( i = 1 ; i < size ; ++i )
	{
		if ( t[k].seVector.IsEqual( t[i].seVector ) )
		{
			t[k].abundance += t[i].abundance ;
			t[i].seVector.Release() ;
		}
		else
		{
			++k ;
			if ( i != k )
				t[k] = t[i] ;
		}
	}
	t.resize( k + 1 ) ;
}

void TranscriptDecider::EnumerateTranscript( int tag, int visit[], int vcnt, struct _subexon *subexons, SubexonCorrelation &correlation, double correlationScore, struct _transcript *alltranscripts, int &atcnt )
{
	int i ;
	visit[ vcnt ] = tag ;
	//printf( "%s: %d %d %d\n", __func__, vcnt, tag, subexons[tag].nextCnt ) ;
	// Compute the correlation score
	double minCor = correlationScore ;
	for ( i = 0 ; i < vcnt - 1 ; ++i )
	{
		double tmp = correlation.Query( visit[i], visit[vcnt] ) ;
		if ( tmp < minCor ) 
			minCor = tmp ;
	}

	if ( subexons[tag].canBeEnd )
	{
		struct _transcript &txpt = alltranscripts[atcnt] ;
		for ( i = 0 ; i <= vcnt ; ++i )
			txpt.seVector.Set( visit[i] ) ;

		txpt.first = visit[0] ;
		txpt.last = visit[vcnt] ;
		txpt.partial = false ;
		txpt.correlationScore = minCor ;
		++atcnt ;
	}

	for ( i = 0 ; i < subexons[tag].nextCnt ; ++i )
		EnumerateTranscript( subexons[tag].next[i], visit, vcnt + 1, subexons, correlation, minCor, alltranscripts, atcnt ) ;
}

// Pick the transcripts from the enumerated transcripts.
void TranscriptDecider::PickTranscripts( struct _transcript *alltranscripts, const int &atcnt, Constraints &constraints, 
		SubexonCorrelation &seCorrelation, std::vector<struct _transcript> &transcripts ) 
{
	int i, j ;
	std::vector<int> chosen ;
	std::vector<struct _matePairConstraint> &tc = constraints.matePairs ;
	int tcCnt = tc.size() ; // transcript constraints
	if ( tcCnt == 0 )
		return ;
	double inf = -1 ; // infinity
	int coalesceThreshold = 1024 ;
	BitTable *btable = new BitTable[ atcnt ] ; 
	for ( i = 0 ; i < atcnt ; ++i )
		btable[i].Init( tcCnt ) ;
	for ( j = 0 ; j < tcCnt ; ++j )
	{
		int a = constraints.matePairs[j].i ;
		int b = constraints.matePairs[j].j ;

		if ( constraints.constraints[a].support > inf )
			inf = constraints.constraints[a].support ;
		if ( constraints.constraints[b].support > inf )
			inf = constraints.constraints[b].support ;

		tc[j].abundance = tc[j].normAbund ;
	}
	++inf ;
	bool btableSet = false ;
	for ( i = 0 ; i < atcnt ; ++i )
	{
		for ( j = 0 ; j < tcCnt ; ++j )
		{
			int a = tc[j].i ;
			int b = tc[j].j ;

			//printf( "try set btble[ %d ].Set( %d ): %d %d\n", i, j, a, b ) ;
			if ( IsConstraintInTranscript( alltranscripts[i], constraints.constraints[a] ) == 1 
					&& IsConstraintInTranscript( alltranscripts[i], constraints.constraints[b] ) == 1 )
			{
				//printf( "set btble[ %d ].Set( %d ): %d %d\n", i, j, a, b ) ;
				btable[i].Set( j ) ;
				btableSet = true ;
			}
		}
	}
	if ( btableSet == false )
	{
		for ( i = 0 ; i < atcnt ; ++i )
			btable[i].Release() ;
		delete[] btable ;
	}

	double maxAbundance = -1 ; // The abundance of the most-abundant transcript
	for ( i = 0 ; i < atcnt ; ++i )
	{
		double value = inf ;
		int tag = -1 ;
		for ( j = 0 ; j < tcCnt ; ++j )
		{
			if ( btable[i].Test(j) && tc[j].abundance > 0 )
			{	
				if ( tc[j].abundance < value )
				{
					value = tc[j].abundance ;
					tag = j ;
				}
			}
		}
		if ( tag == -1 ) 
			value = 1e-6 ;
		if ( value > maxAbundance )
			maxAbundance = value ;
	}

	while ( 1 )
	{
		double max = -1 ;
		int maxtag = -1 ;
		int maxcnt = -1 ;

		// Find the optimal candidate.
		for ( i = 0 ; i < atcnt ; ++i )
		{
			double value = inf ;
			double cnt = 0 ;
			int allCnt = 0 ;

			for ( j = 0 ; j < tcCnt ; ++j )
			{
				if ( btable[i].Test( j ) )
				{
					if ( tc[j].normAbund < value )
						value = tc[j].normAbund ;

					if ( tc[j].abundance > 0 )
					{		
						cnt += tc[j].effectiveCount ;
					}
				}
			}

			if ( value == inf )
				value = 1e-6 ;

			double score = ComputeScore( cnt, value, maxAbundance, alltranscripts[i].correlationScore ) ;
			if ( cnt > maxcnt )
				maxcnt = cnt ;

			if ( score > max )
			{
				max = score ;
				maxtag = i ;
			}
		}

		if ( maxcnt == 0 || maxtag == -1 )
			break ;

		// Update the abundance for each transcript.	
		double update = inf ;
		for ( j = 0 ; j < tcCnt ; ++j )
		{
			if ( btable[ maxtag ].Test( j ) && tc[j].abundance > 0 && 
					tc[j].abundance < update )
			{
				update = tc[j].abundance ;	
			}
		}

		struct _transcript nt ;
		//printf( "maxtag=%d %lf\n", maxtag, update ) ;
		nt.seVector.Duplicate( alltranscripts[ maxtag ].seVector ) ;
		nt.first = alltranscripts[maxtag].first ;
		nt.last = alltranscripts[maxtag].last ;
		nt.abundance = 0 ; 
		for ( j = 0 ; j < tcCnt ; ++j )
		{
			if ( btable[maxtag].Test( j ) && tc[j].abundance > 0 )
			{
				tc[j].abundance -= 1 * update ;
				double factor = tc[j].effectiveCount ;

				nt.abundance += ( tc[j].support * update / tc[j].normAbund * factor ) ;
			}

			if ( tc[j].abundance < 0 )
				tc[j].abundance = 0 ;
		}

		transcripts.push_back( nt ) ;
		if ( transcripts.size() >= transcripts.capacity() && (int)transcripts.size() >= coalesceThreshold )
		{
			CoalesceSameTranscripts( transcripts ) ;
			if ( transcripts.size() >= transcripts.capacity() / 2 )
				coalesceThreshold *= 2 ;
		}
	}

	// Clean up the transcripts.
	CoalesceSameTranscripts( transcripts ) ;

	// Release the memory of btable.
	for ( i = 0 ; i < atcnt ; ++i )
		btable[i].Release() ;
	delete[] btable ;
}

int TranscriptDecider::Solve( struct _subexon *subexons, int seCnt, std::vector<Constraints> &constraints, SubexonCorrelation &subexonCorrelation )
{
	int i, j ;
	int cnt = 0 ;
	int *f = new int[seCnt] ; // this is a general buffer for a type of usage.	
	bool useDP = false ;

	compatibleTestVectorT.Init( seCnt ) ; // this is the bittable used in compatible test function.	
	compatibleTestVectorC.Init( seCnt ) ;

	for ( i = 0 ; i < seCnt ; ++i )
	{
		subexons[i].canBeStart = subexons[i].canBeEnd = false ;

		if ( subexons[i].prevCnt == 0 )
			subexons[i].canBeStart = true ;
		else if ( subexons[i].leftClassifier < canBeSoftBoundaryThreshold )
			subexons[i].canBeStart = true ;

		if ( subexons[i].nextCnt == 0 )
			subexons[i].canBeEnd = true ;
		else if ( subexons[i].rightClassifier < canBeSoftBoundaryThreshold )
			subexons[i].canBeEnd = true ;
	}

	geneId = new int[ seCnt ] ;
	memset( geneId, -1, sizeof( int ) * seCnt ) ;
	int baseGeneId = usedGeneId ;
	for ( i = 0 ; i < seCnt ; ++i )
	{
		if ( geneId[i] == -1 )
		{
			SetGeneId( i, subexons, usedGeneId ) ;
			++usedGeneId ;
		}
	}

	cnt = 0 ;
	memset( f, -1, sizeof( int ) * seCnt ) ;
	for ( i = 0 ; i < seCnt ; ++i )
	{
		if ( subexons[i].canBeStart )	
		{
			cnt += SubTranscriptCount( i, subexons, f ) ;
		}
	}
	if ( cnt <= USE_DP )
	{
		for ( i = 0 ; i < seCnt ; ++i )
			if ( f[i] > USE_DP )
			{
				useDP = true ;
				break ;
			}
	}
	else
		useDP = true ;
	if ( !useDP )
	{
		for ( i = 0 ; i < sampleCnt ; ++i )
			if ( constraints[i].matePairs.size() > 3000 && cnt > 200 )
			{
				useDP = true ;
				break ;
			}
	}

	int atCnt = cnt ;

	if ( !useDP )
	{
		struct _transcript *alltranscripts = new struct _transcript[atCnt] ;

		for ( i = 0 ; i < atCnt ; ++i )
		{
			alltranscripts[i].seVector.Init( seCnt ) ; 
			alltranscripts[i].correlationScore = 1 ;
		}

		atCnt = 0 ;
		for ( i = 0 ; i < seCnt ; ++i )
		{
			if ( subexons[i].canBeStart )
				EnumerateTranscript( i, f, 0, subexons, subexonCorrelation, 1, alltranscripts, atCnt ) ;
		}
		//printf( "transcript cnt: %d\n", atCnt ) ;
		//printf( "%d %d\n", alltranscripts[0].seVector.Test( 1 ), constraints[0].matePairs.size() ) ;
		transcriptId = new int[usedGeneId - baseGeneId] ;
		for ( i = 0 ; i < sampleCnt ; ++i )
		{
			std::vector<struct _transcript> predTranscripts ;
			PickTranscripts( alltranscripts, atCnt, constraints[i], subexonCorrelation, predTranscripts ) ;
			int size = predTranscripts.size() ;
			InitTranscriptId( baseGeneId, usedGeneId ) ;
			for ( j = 0 ; j < size ; ++j )
			{
				OutputTranscript( outputFPs[i], baseGeneId, subexons, predTranscripts[j] ) ;
			}
			for ( j = 0 ; j < size ; ++j )
				predTranscripts[j].seVector.Release() ;
		}
		delete []transcriptId ;

		for ( i = 0 ; i < atCnt ; ++i )
			alltranscripts[i].seVector.Release() ;
		delete[] alltranscripts ;
	}
	else
	{
		;
	}

	delete[] geneId ;
	return 0 ;	
}

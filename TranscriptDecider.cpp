#include "TranscriptDecider.hpp"

#define USE_DP 200000
#define HASH_MAX 1000003 


void TranscriptDecider::OutputTranscript( int sampleId, struct _subexon *subexons, struct _transcript &transcript )
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
			/*int nextCnt = subexons[ subexonInd[i] ].nextCnt ;
			if ( nextCnt == 0 )
				continue ;

			for ( j = 0 ; j < nextCnt ; ++j )
			{
				int a = subexons[ subexonInd[i] ].next[j] ;
				if ( subexonInd[i + 1] == a 
					&& subexons[ subexonInd[i] ].end + 1 < subexons[a].start ) // avoid the case like ..(...[...
				{
					break ;
				}
			}
			if ( j < nextCnt )*/

			if ( subexons[ subexonInd[i] ].end + 1 < subexons[ subexonInd[i + 1] ].start )
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
	char prefix[10] = "" ;
	struct _subexon *catSubexons = new struct _subexon[ size + 1 ] ;
	// Concatenate adjacent subexons 
	catSubexons[0] = subexons[ subexonInd[0] ] ;
	j = 1 ;
	for ( i = 1 ; i < size ; ++i )
	{
		if ( subexons[ subexonInd[i] ].start == catSubexons[j - 1].end + 1 )
		{
			catSubexons[j - 1].end = subexons[ subexonInd[i] ].end ;
		}
		else		
		{
			catSubexons[j] = subexons[ subexonInd[i] ] ;
			++j ;
		}
	}
	size = j ;
	
	int gid = GetTranscriptGeneId( subexonInd, subexons ) ;
	if ( numThreads <= 1 )
	{
		fprintf( outputFPs[sampleId], "%s\tCLASSES\ttranscript\t%d\t%d\t1000\t%s\t.\tgene_id \"%s%s.%d\"; transcript_id \"%s%s.%d.%d\"; Abundance \"%.6lf\";\n",
				chrom, catSubexons[0].start + 1, catSubexons[size - 1].end + 1, strand,
				prefix, chrom, gid,
				prefix, chrom, gid, transcriptId[ gid - baseGeneId ], transcript.FPKM ) ;
		for ( i = 0 ; i < size ; ++i )
		{
			fprintf( outputFPs[ sampleId ], "%s\tCLASSES\texon\t%d\t%d\t1000\t%s\t.\tgene_id \"%s%s.%d\"; "
					"transcript_id \"%s%s.%d.%d\"; exon_number \"%d\"; Abundance \"%.6lf\"\n",
					chrom, catSubexons[i].start + 1, catSubexons[i].end + 1, strand,
					prefix, chrom, gid,
					prefix, chrom, gid, transcriptId[ gid - baseGeneId ],
					i + 1, transcript.FPKM ) ;
		}
	}
	else
	{
		struct _outputTranscript t ;	
		t.chrId = subexons[0].chrId ;
		t.geneId = gid ;
		t.transcriptId = transcriptId[ gid - baseGeneId ] ;
		t.FPKM = transcript.FPKM ;
		t.sampleId = sampleId ;
		t.exons = new struct _pair32[size] ;
		for ( i = 0 ; i < size ; ++i )
		{
			t.exons[i].a = catSubexons[i].start + 1 ;
			t.exons[i].b = catSubexons[i].end + 1 ;
		}
		t.ecnt = size ;
		t.strand = strand[0] ;

		outputHandler->Add( t ) ;
	}
	++transcriptId[ gid - baseGeneId ] ;

	delete[] catSubexons ;
}


int TranscriptDecider::GetTranscriptGeneId( std::vector<int> &subexonInd, struct _subexon *subexons )
{
	int i ;
	int size = subexonInd.size() ;

	for ( i = 0 ; i < size ; ++i )
		if ( subexons[ subexonInd[i] ].geneId != -2  )
			return  subexons[ subexonInd[i] ].geneId ;
	
	// Some extreme case, where all the regions are mixture regions.
	for ( i = 0 ; i < size - 1 ; ++i )
		if ( subexons[ subexonInd[i] ].end + 1 < subexons[ subexonInd[i + 1] ].start )
		{
			return defaultGeneId[ ( subexons[ subexonInd[i] ].rightStrand + 1 ) / 2 ] ;
		}
	return defaultGeneId[0] ;
}

int TranscriptDecider::GetTranscriptGeneId( struct _transcript &t, struct _subexon *subexons )
{
	if ( subexons[ t.first ].geneId != -2 )
		return subexons[ t.first ].geneId ;
	if ( subexons[ t.last ].geneId != -2 )
		return subexons[ t.last ].geneId ;
	std::vector<int> subexonInd ;
	t.seVector.GetOnesIndices( subexonInd ) ;
	return GetTranscriptGeneId( subexonInd, subexons ) ;
}

void TranscriptDecider::InitTranscriptId()
{
	int i ;
	for ( i = 0 ; i < usedGeneId - baseGeneId ; ++i )
		transcriptId[i] = 0 ;
}

bool TranscriptDecider::IsStartOfMixtureStrandRegion( int tag, struct _subexon *subexons, int seCnt ) 
{
	int j, k ;
	int leftStrandCnt[2] = {0, 0}, rightStrandCnt[2] = {0, 0};
	for ( j = tag + 1 ; j < seCnt ; ++j )
		if ( subexons[j].start > subexons[j - 1].end + 1 )
			break ;
	
	for ( k = tag ; k < j ; ++k )
	{
		if ( subexons[k].leftStrand != 0 )
			++leftStrandCnt[ ( subexons[k].leftStrand + 1 ) / 2 ] ;
		if ( subexons[k].rightStrand != 0 )
			++rightStrandCnt[ ( subexons[k].rightStrand + 1 ) / 2 ] ;
	}

	if ( rightStrandCnt[0] > 0 && leftStrandCnt[0] == 0 && leftStrandCnt[1] > 0 )
		return true ;
	if ( rightStrandCnt[1] > 0 && leftStrandCnt[1] == 0 && leftStrandCnt[0] > 0 )
		return true ;
	return false ;
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
	if ( c.last > transcript.last )
	{
		compatibleTestVectorC.MaskRegionOutside( s, e ) ;
	}
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

int TranscriptDecider::IsConstraintInTranscriptDebug( struct _transcript transcript, struct _constraint &c ) 
{
	//printf( "%d %d, %d %d\n", c.first, c.last, transcript.first, transcript.last ) ;
	if ( c.first < transcript.first || c.first > transcript.last ) // no overlap or starts too early.
		return 0 ; 
	printf( "hi\n" ) ;
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
	compatibleTestVectorT.Print() ;
	compatibleTestVectorC.Print() ;

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
	if ( t.size() == 0 )
		return ;

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

void TranscriptDecider::EnumerateTranscript( int tag, int strand, int visit[], int vcnt, struct _subexon *subexons, SubexonCorrelation &correlation, double correlationScore, std::vector<struct _transcript> &alltranscripts, int &atcnt )
{
	int i ;
	visit[ vcnt ] = tag ;
	//printf( "%s: %d %d %d %d. %d %d\n", __func__, vcnt, tag, subexons[tag].nextCnt, strand, subexons[tag].start, subexons[tag].end ) ;
	// Compute the correlation score
	double minCor = correlationScore ;
	for ( i = 0 ; i < vcnt ; ++i )
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
		
		//printf( "%lf %d %d ", txpt.correlationScore, vcnt, visit[0] ) ;
		//txpt.seVector.Print() ;
		++atcnt ;
	}

	for ( i = 0 ; i < subexons[tag].nextCnt ; ++i )
	{
		int a = subexons[tag].next[i] ;
		if ( !SubexonGraph::IsSameStrand( subexons[tag].rightStrand, strand ) 
			&& subexons[a].start > subexons[tag].end + 1 )
			continue ;
		int backupStrand = strand ;
		if ( subexons[a].start > subexons[tag].end + 1 && strand == 0 )
			strand = subexons[tag].rightStrand ;
		EnumerateTranscript( subexons[tag].next[i], strand, visit, vcnt + 1, subexons, correlation, minCor, alltranscripts, atcnt ) ;
		strand = backupStrand ;
	}
}

void TranscriptDecider::SearchSubTranscript( int tag, int strand, int parents[], int pcnt, struct _dp &pdp, int visit[], int vcnt, int extends[], int extendCnt,
std::vector<struct _constraint> &tc, int tcStartInd, struct _dpAttribute &attr ) 
{
	int i ;
	int size ;
	double cover ;
	bool keepSearch = true ;
	bool belowMin = false ;
	
	struct _subexon *subexons = attr.subexons ;

	visit[vcnt] = tag ;
	++vcnt ;
	struct _dp visitdp ;
	
	visitdp.cover = -1 ;
	
	struct _transcript &subTxpt = attr.bufferTxpt ;
	subTxpt.seVector.Reset() ;
	for ( i = 0 ; i < pcnt ; ++i ) 
		subTxpt.seVector.Set( parents[i] ) ;
	subTxpt.first = parents[0] ;
	subTxpt.last = parents[ pcnt - 1] ;
	for ( i = 0 ; i < vcnt ; ++i )
		subTxpt.seVector.Set( visit[i] ) ;
	subTxpt.partial = true ;

	// Adjust the extendsCnt
	for ( i = extendCnt - 1 ; i >= 0 ; --i )
		if ( IsConstraintInTranscript( subTxpt, tc[ extends[i] ] ) != 0 )
			break ;
	extendCnt = i + 1 ;
	
	// If the extension ends.
	subTxpt.partial = false ;
	if ( subexons[tag].nextCnt > 0 && ( extendCnt == 0 || tag >= tc[ extends[ extendCnt - 1 ] ].last ) )
	{
		// Solve the subtranscript beginning with visit.
		// Now we got the optimal transcript for visit. 
		visitdp = SolveSubTranscript( visit, vcnt, strand, tc, tcStartInd, attr ) ;	
		keepSearch = false ;
	}
	//printf( "%s %d %d: visitdp.cover=%lf\n", __func__, parents[0], tag, visitdp.cover ) ;

	// the constraints across the parents and visit.
	size = tc.size() ;
	if ( visitdp.cover >= 0 )
	{
		cover = visitdp.cover ;
		// Reset the subTxpt, since its content is modofitied in SolveSubTxpt called above.
		subTxpt.seVector.Reset() ;
		for ( i = 0 ; i < pcnt ; ++i ) 
			subTxpt.seVector.Set( parents[i] ) ;
		subTxpt.seVector.Or( visitdp.seVector ) ;
		subTxpt.first = parents[0] ;
		subTxpt.last = visitdp.last ;
		subTxpt.partial = false ;

		for ( i = tcStartInd ; i < size ; ++i )		
		{
			if ( tc[i].first > parents[ pcnt - 1] )
				break ;

			if ( IsConstraintInTranscript( subTxpt, tc[i] ) == 1 )
			{
				if ( tc[i].normAbund <= attr.minAbundance )
				{
					belowMin = true ;
					cover = -2 ;
					break ;
				}

				if ( tc[i].abundance <= 0 )
					continue ;

				if ( attr.forAbundance )
				{
					if ( tc[i].normAbund < cover || cover == 0 )
						cover = tc[i].normAbund ;
				}
				else
				{
					++cover ;
				}
			}
		}
		if ( belowMin && pdp.cover == -1 )	
		{
			pdp.cover = -2 ;
			pdp.strand = strand ;
		}
		else if ( cover > pdp.cover )
		{
			pdp.cover = cover ;
			pdp.seVector.Assign( subTxpt.seVector ) ;
			pdp.first = subTxpt.first ;
			pdp.last = subTxpt.last ;
			pdp.strand = strand ;
		}
	}
	else if ( visitdp.cover == -2 && pdp.cover == -1 )
	{
		pdp.cover = -2 ;
		pdp.strand = strand ;
	}

	if ( subexons[tag].canBeEnd && ( visitdp.cover < 0 || attr.forAbundance ) ) 
	// This works is because that the extension always covers more constraints. So we only go this branch if the extension does not work
	// and it goes this branch if it violates minAbundance
	// But we need to go here when we want to compute the maxAbundance transcript.
	// This part also works as the exit point of the recurive function.
	{
		bool belowMin = false ;
		subTxpt.seVector.Reset() ;
		for ( i = 0 ; i < pcnt ; ++i ) 
			subTxpt.seVector.Set( parents[i] ) ;
		for ( i = 0 ; i < vcnt ; ++i )
			subTxpt.seVector.Set( visit[i] ) ;
		subTxpt.first = parents[0] ;
		subTxpt.last = visit[ vcnt - 1] ;
		subTxpt.partial = false ;

		cover = 0 ;
		for ( i = tcStartInd ; i < size ; ++i )		
		{ 
			// note that the value is parents[ pcnt - 1], because  
			// in above the part of "visit" is computed in SolveSubTranscript( visit ).
			if ( tc[i].first > visit[ vcnt - 1] )
				break ;
			if ( IsConstraintInTranscript( subTxpt, tc[i] ) == 1 )
			{
				if ( tc[i].normAbund <= attr.minAbundance )
				{
					belowMin = true ;
					cover = -2 ;
					break ;
				}

				if ( tc[i].abundance <= 0 )
					continue ;
				if ( attr.forAbundance )
				{
					if ( tc[i].normAbund < cover || cover == 0 )	
						cover = tc[i].normAbund ;
				}
				else
				{
					++cover ;
				}
			}
		}

		if ( belowMin && pdp.cover == -1 )	
		{
			pdp.cover = -2 ;
			pdp.strand = strand ;
		}
		else if ( cover > pdp.cover )
		{
			pdp.cover = cover ;
			pdp.seVector.Assign( subTxpt.seVector ) ;
			pdp.first = subTxpt.first ;
			pdp.last = subTxpt.last ;
			pdp.strand = strand ;
		}
	}
	//printf( "%s %d: pdp.cover=%lf\n", __func__, tag, pdp.cover ) ;

	// keep searching.
	if ( keepSearch )	
	{
		for ( i = 0 ; i < subexons[tag].nextCnt ; ++i )
		{
			int b = subexons[tag].next[i] ;
			if ( SubexonGraph::IsSameStrand( subexons[tag].rightStrand, strand ) ||
					subexons[b].start == subexons[tag].end + 1 )		
			{
				int backupStrand = strand ;
				if ( subexons[b].start > subexons[tag].end + 1 ) 
					strand = subexons[tag].rightStrand ;

				SearchSubTranscript( subexons[tag].next[i], strand, parents, pcnt, pdp, visit, vcnt, 
						extends, extendCnt, tc, tcStartInd, attr ) ;
				strand = backupStrand ;
			}
		}

	}

	return ;
}

struct _dp TranscriptDecider::SolveSubTranscript( int visit[], int vcnt, int strand, std::vector<struct _constraint> &tc, int tcStartInd, struct _dpAttribute &attr ) 
{
	int i ;
	int size ;
	/*printf( "%s: ", __func__ ) ;	
	for ( i = 0 ; i < vcnt ; ++i )
		printf( "%d ", visit[i] ) ;
	printf( ": %lf %d %d", attr.f1[ visit[0] ].cover, attr.f1[ visit[0] ].timeStamp, attr.timeStamp ) ;
	printf( "\n" ) ;*/
	// Test whether it is stored in dp 
	if ( vcnt == 1 )
	{
		if ( attr.f1[ visit[0] ].cover != -1 && attr.f1[ visit[0] ].strand == strand && ( attr.f1[ visit[0] ].timeStamp == attr.timeStamp  || 
			( attr.f1[ visit[0] ].minAbundance < attr.minAbundance && attr.f1[visit[0]].cover == -2 ) ) )
			return attr.f1[ visit[0] ] ;
	}
	else if ( vcnt == 2 )
	{
		int a = visit[0] ;
		int b = visit[1] ;
		
		if ( attr.f2[a][b].cover != -2 && attr.f2[a][b].strand == strand && ( attr.f2[a][b].timeStamp == attr.timeStamp || 
			( attr.f2[a][b].minAbundance < attr.minAbundance && attr.f2[a][b].cover == -2 ) ) )
		{
			return attr.f2[a][b] ;
		}
	}
	else
	{
		int key = 0 ;	
		for ( i = 0 ; i < vcnt ; ++i )
			key = ( key * 17 + visit[i] ) % HASH_MAX ;

		if ( attr.hash[key].cover != -1 && attr.hash[key].cnt == vcnt && attr.hash[key].strand == strand && 
			( attr.hash[key].timeStamp == attr.timeStamp || 
				( attr.hash[key].minAbundance < attr.minAbundance && attr.hash[key].cover == -2 ) ) )
		{
			struct _transcript subTxpt = attr.bufferTxpt ;
			subTxpt.seVector.Reset() ;
			for ( i = 0 ; i < vcnt ; ++i )
				subTxpt.seVector.Set( visit[i] ) ;
			subTxpt.seVector.Xor( attr.hash[key].seVector ) ;
			subTxpt.seVector.MaskRegionOutside( visit[0], visit[ vcnt - 1] ) ;
			if ( subTxpt.seVector.IsAllZero() )
			{
				return attr.hash[key] ;
			}
		}
	}
	// adjust tcStartInd 
	size = tc.size() ;
	for ( i = tcStartInd ; i < size ; ++i )
		if ( tc[i].first >= visit[0] )
			break ;
	tcStartInd = i ;


	struct _subexon *subexons = attr.subexons ;
	struct _dp visitdp ;
	visitdp.seVector.Init( attr.seCnt ) ;
	visitdp.cover = -1 ;

	struct _transcript &subTxpt = attr.bufferTxpt ;
	// This happens when it is called from PickTranscriptsByDP, the first subexon might be the end.
	subTxpt.seVector.Reset() ;	
	for ( i = 0 ; i < vcnt ; ++i )
		subTxpt.seVector.Set( visit[i] ) ;
	subTxpt.first = visit[0] ;
	subTxpt.last = visit[vcnt - 1] ;

	if ( subexons[ visit[vcnt - 1] ].canBeEnd )
	{
		subTxpt.partial = false ;
		double cover = 0 ;
		for ( i = tcStartInd ; i < size ; ++i )		
		{
			if ( tc[i].first > subTxpt.last )
				break ;

			if ( IsConstraintInTranscript( subTxpt, tc[i] ) == 1 )
			{
				if ( tc[i].normAbund <= attr.minAbundance )
				{
					cover = -2 ;
					break ;
				}

				if ( tc[i].abundance <= 0 )
					continue ;
				if ( attr.forAbundance )
				{
					if ( tc[i].normAbund < cover || cover == 0 )	
						cover = tc[i].normAbund ;
				}
				else
					++cover ;	
			}
		}

		visitdp.seVector.Assign( subTxpt.seVector ) ;		
		visitdp.cover = cover ;
		visitdp.first = subTxpt.first ;
		visitdp.last = subTxpt.last ;
		visitdp.strand = strand ;
	}
	
	// Now we extend.
	size = tc.size() ;
	int *extends = new int[tc.size() - tcStartInd + 1] ;
	int extendCnt = 0 ;
	subTxpt.partial = true ;
	for ( i = tcStartInd ; i < size ; ++i )
	{
		if ( tc[i].first > subTxpt.last )
			break ;
		if ( IsConstraintInTranscript( subTxpt, tc[i] ) == 2 )
		{
			extends[extendCnt] = i ;
			++extendCnt ;
		}
	}

	// Sort the extend by the index of the last subexon.
	if ( extendCnt > 0 )
	{
		struct _pair32 *extendsPairs = new struct _pair32[extendCnt] ;

		for ( i = 0 ; i < extendCnt ; ++i )
		{
			extendsPairs[i].a = extends[i] ;
			extendsPairs[i].b = tc[ extends[i] ].last ;
		}
		qsort( extendsPairs, extendCnt, sizeof( struct _pair32 ), CompExtendsPairs ) ;

		for ( i = 0 ; i < extendCnt ; ++i )
			extends[i] = extendsPairs[i].a ;

		delete[] extendsPairs ;
	}

	size = subexons[ visit[vcnt - 1] ].nextCnt ;
	int nextvCnt = 1 ;
	if ( extendCnt > 0 && tc[ extends[ extendCnt - 1 ] ].last - visit[ vcnt - 1 ] > 1 )
		nextvCnt = tc[ extends[ extendCnt - 1 ] ].last - visit[ vcnt - 1 ] ;
	int *nextv = new int[ nextvCnt ] ;
	for ( i = 0 ; i < size ; ++i )
	{
		int a = visit[vcnt - 1] ;
		int b = subexons[a].next[i] ;
		if ( SubexonGraph::IsSameStrand( subexons[a].rightStrand, strand ) ||
			subexons[b].start == subexons[a].end + 1 )		
		{
			int backupStrand = strand ;
			if ( subexons[b].start > subexons[a].end + 1 ) 
				strand = subexons[a].rightStrand ;
			SearchSubTranscript( subexons[ visit[vcnt - 1] ].next[i], strand, visit, vcnt, visitdp, nextv, 0, extends, extendCnt, tc, tcStartInd, attr ) ;		
			strand = backupStrand ;

		}
	}
	//printf( "%s %d(%d) %d %d %d: %lf\n", __func__, visit[0], subexons[ visit[vcnt - 1] ].canBeEnd, size, extendCnt, strand, visitdp.cover ) ;	
	delete[] nextv ;
	delete[] extends ;

	// store the result in the dp structure.
	// We return the structure stored in dp to simplify the memory access pattern.
	// In other words, we assume the structure returned from this function always uses the memory from attr.dp 
	if ( vcnt == 1 )
	{
		SetDpContent( attr.f1[ visit[0] ], visitdp, attr ) ;
		visitdp.seVector.Release() ;
		return attr.f1[ visit[0] ] ;
	}
	else if ( vcnt == 2 )
	{
		SetDpContent( attr.f2[ visit[0] ][ visit[1] ], visitdp, attr ) ;
		visitdp.seVector.Release() ;
		return attr.f2[ visit[0] ][ visit[1] ] ;
	}
	else
	{
		int key = 0 ;	
		for ( i = 0 ; i < vcnt ; ++i )
			key = ( key * 17 + visit[i] ) % HASH_MAX ;
		SetDpContent( attr.hash[key], visitdp, attr ) ;	
		attr.hash[key].cnt = vcnt ;
		visitdp.seVector.Release() ;
		return attr.hash[key] ;
	}
}

void TranscriptDecider::PickTranscriptsByDP( struct _subexon *subexons, int seCnt, Constraints &constraints, SubexonCorrelation &correlation, struct _dpAttribute &attr, std::vector<struct _transcript> &alltranscripts )
{
	int i, j, k ;
	
	std::vector<struct _transcript> transcripts ;
	std::vector<struct _constraint> &tc = constraints.constraints ;
	int tcCnt = tc.size() ;
	int coalesceThreshold = 1024 ;

	printf( "tcCnt=%d\n", tcCnt ) ;

	attr.timeStamp = 1 ;
	attr.bufferTxpt.seVector.Init( seCnt ) ;
	attr.subexons = subexons ;
	attr.seCnt = seCnt ;

	double maxAbundance = -1 ;
	// Initialize the dp data structure
	/*memset( attr.f1, -1, sizeof( struct _dp ) * seCnt ) ;
	for ( i = 0 ; i < seCnt ; ++i )
		memset( attr.f2[i], -1, sizeof( struct _dp ) * seCnt ) ;
	memset( attr.hash, -1, sizeof( struct _dp ) * HASH_MAX ) ;*/
	for ( i = 0 ; i < seCnt ; ++i )
		ResetDpContent( attr.f1[i] ) ;
	for ( i = 0 ; i < seCnt ; ++i )
		for ( j = i ; j < seCnt ; ++j )
			ResetDpContent( attr.f2[i][j] ) ;
	for ( i = 0 ; i < HASH_MAX ; ++i )
		ResetDpContent( attr.hash[i] ) ;

	// Find the max abundance 
	attr.forAbundance = true ;
	attr.minAbundance = 0 ;
	for ( i = 0 ; i < seCnt ; ++i )
	{
		if ( subexons[i].canBeStart )	
		{
			int visit[1] = {i} ;
			struct _dp tmp ;
			
			tmp = SolveSubTranscript( visit, 1, 0, tc, 0, attr ) ;
			
			if ( tmp.cover > maxAbundance )
				maxAbundance = tmp.cover ;
		}
	}
	//printf( "maxAbundance=%lf\n", maxAbundance ) ;
	
	// Pick the transcripts
	// Notice that by the logic in SearchSubTxpt and SolveSubTxpt, we don't need to reinitialize the data structure.
	attr.forAbundance = false ;
	int *coveredTc = new int[tcCnt] ;
	int coveredTcCnt ;
	struct _dp maxCoverDp ;
	struct _dp bestDp ;

	maxCoverDp.seVector.Init( seCnt ) ;
	bestDp.seVector.Init( seCnt ) ;
	while ( 1 )
	{
		double bestScore ;
		
		// iterately assign constraints
		attr.minAbundance = 0 ;	
		
		// Find the best candidate transcript.
		bestDp.cover = -1 ;
		bestScore = -1 ;
		while ( 1 )
		{
			// iterate the change of minAbundance
			maxCoverDp.cover = -1 ;
			++attr.timeStamp ;
			for ( i = 0 ; i < seCnt ; ++i )		
			{
				if ( subexons[i].canBeStart == false )
					continue ;
				int visit[1] = {i} ;
				struct _dp tmp ;
				tmp = SolveSubTranscript( visit, 1, 0, tc, 0, attr ) ;

				if ( tmp.cover > maxCoverDp.cover && tmp.cover > 0 )
				{
					SetDpContent( maxCoverDp, tmp, attr ) ;
				}
				//if ( subexons[i].start == 6870264 || subexons[i].start == 6872237 )
				//	printf( "%d: %lf\n", i, tmp.cover ) ;
			}
			
			if ( maxCoverDp.cover == -1 )
				break ;
			// the abundance for the max cover txpt.
			double min = -1 ;
			struct _transcript &subTxpt = attr.bufferTxpt ;
			subTxpt.seVector.Assign( maxCoverDp.seVector ) ;
			subTxpt.first = maxCoverDp.first ;
			subTxpt.last = maxCoverDp.last ;

			for ( i = 0 ; i < tcCnt ; ++i )
			{
				if ( IsConstraintInTranscript( subTxpt, tc[i] ) == 1 )	
				{
					if ( tc[i].normAbund < min || min == -1 )
						min = tc[i].normAbund ;
				}
			}
			
			double score = ComputeScore( maxCoverDp.cover, 1.0, min, maxAbundance, 0 ) ;
			if ( bestScore == -1 || score > bestScore )	
			{
				bestScore = score ;
				SetDpContent( bestDp, maxCoverDp, attr ) ;
			}
			//printf( "normAbund=%lf maxCoverDp.cover=%lf score=%lf\n", min, maxCoverDp.cover, score ) ;
			attr.minAbundance = min ;
		}

		if ( bestDp.cover == -1 )
			break ;
		// Assign the constraints.
		coveredTcCnt = 0 ;
		double update = -1 ;
		struct _transcript &subTxpt = attr.bufferTxpt ;
		subTxpt.seVector.Assign( bestDp.seVector ) ;
		subTxpt.first = bestDp.first ;
		subTxpt.last = bestDp.last ;
		subTxpt.partial = false ;
		for ( i = 0 ; i < tcCnt ; ++i )
		{
			if ( IsConstraintInTranscript( subTxpt, tc[i] ) == 1 )
			{
				if ( tc[i].abundance > 0 && 
					( tc[i].abundance < update || update == -1 ) )
				{
					update = tc[i].abundance ;		
				}
				coveredTc[ coveredTcCnt ] = i ;
				++coveredTcCnt ;
			}
			/*else
			{
				printf( "%d: ", i ) ;
				tc[i].vector.Print() ;
				if ( i == 127 )
				{
					printf( "begin debug:\n" ) ;
					IsConstraintInTranscriptDebug( subTxpt, tc[i] ) ;
				}
			}*/
		}

		/*printf( "update=%lf %d %d. %d %d %d\n", update, coveredTcCnt, tcCnt, 
				bestDp.first, bestDp.last, subexons[ bestDp.first ].start ) ;
		bestDp.seVector.Print() ;*/
		
		struct _transcript nt ;
		nt.seVector.Duplicate( bestDp.seVector ) ; 
		nt.first = bestDp.first ;
		nt.last = bestDp.last ;
		nt.abundance = 0 ;
		for ( i = 0 ; i < coveredTcCnt ; ++i )
		{
			j = coveredTc[i] ;
			if ( tc[j].abundance > 0 )	
			{
				tc[j].abundance -= 1 * update ;
				double factor = 1 ;

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
	CoalesceSameTranscripts( transcripts ) ;
	int size = transcripts.size() ;
	// Compute the correlation score
	for ( i = 0 ; i < size ; ++i )
	{
		std::vector<int> subexonInd ;
		transcripts[i].seVector.GetOnesIndices( subexonInd ) ;
		double cor = 2.0 ;
		int s = subexonInd.size() ;
		for ( j = 0 ; j < s ; ++j )
			for ( k = j + 1 ; j < s ; ++j )
			{
				double tmp = correlation.Query( subexonInd[j], subexonInd[k] ) ;
				if ( tmp < cor )					
					cor = tmp ;
			}
		if ( cor > 1 )
			cor = 0 ;
		transcripts[i].correlationScore = cor ;
	}

	// store the result
	for ( i = 0 ; i < size ; ++i )
		alltranscripts.push_back( transcripts[i] ) ;

	// Release the memory
	attr.bufferTxpt.seVector.Release() ;

	delete[] coveredTc ;	
	maxCoverDp.seVector.Release() ;
	bestDp.seVector.Release() ;
}


// Pick the transcripts from given transcripts.
void TranscriptDecider::PickTranscripts( std::vector<struct _transcript> &alltranscripts, Constraints &constraints, 
		SubexonCorrelation &seCorrelation, std::vector<struct _transcript> &transcripts ) 
{
	int i, j, k ;
	std::vector<int> chosen ;
	std::vector<struct _matePairConstraint> &tc = constraints.matePairs ;
	int atcnt = alltranscripts.size() ;
	int tcCnt = tc.size() ; // transcript constraints
	if ( tcCnt == 0 )
		return ;
	double inf = -1 ; // infinity
	int coalesceThreshold = 1024 ;
	int *transcriptSeCnt = new int[ atcnt ] ;
	double *transcriptAbundance = new double[atcnt] ; // the roughly estimated abundance based on constraints.
	double *avgTranscriptAbundance = new double[atcnt] ; // the average normAbund from the compatible constraints.

	BitTable *btable = new BitTable[ atcnt ] ; 
	//BitTable lowCovSubexon ; // force the abundance to 0 for the transcript contains the subexon.
	double *coveredPortion = new double[atcnt] ;

	memset( avgTranscriptAbundance, 0 ,sizeof( double ) * atcnt ) ;
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
		//printf( "correlation %d: %lf\n", i, alltranscripts[i].correlationScore ) ;
		/*for ( int l = 0 ; l < subexonInd.size() ; ++l )
		{
			for ( int m = l ; m < subexonInd.size() ; ++m )
				printf( "%lf ", seCorrelation.Query( l, m ) ) ;
			printf( "\n" ) ;
		}*/

		for ( j = 0 ; j < tcCnt ; ++j )
		{
			int a = tc[j].i ;
			int b = tc[j].j ;

			//printf( "try set btble[ %d ].Set( %d ): %d %d\n", i, j, a, b ) ;
			//alltranscripts[i].seVector.Print() ;
			//constraints.constraints[a].vector.Print() ;
			//constraints.constraints[b].vector.Print() ;
			if ( IsConstraintInTranscript( alltranscripts[i], constraints.constraints[a] ) == 1 
					&& IsConstraintInTranscript( alltranscripts[i], constraints.constraints[b] ) == 1 )
			{
				//printf( "set btble[ %d ].Set( %d ): %d %d\n", i, j, a, b ) ;
				btable[i].Set( j ) ;
				btableSet = true ;
			}
		}
		transcriptSeCnt[i] = alltranscripts[i].seVector.Count() ;
	}
	if ( btableSet == false )
	{
		for ( i = 0 ; i < atcnt ; ++i )
			btable[i].Release() ;
		delete[] btable ;
		return ;
	}

	double maxAbundance = -1 ; // The abundance of the most-abundant transcript
	if ( atcnt > 0 && alltranscripts[0].abundance == -1 )
	{
		int seCnt = 0 ;
		if ( atcnt > 0 )
			seCnt = alltranscripts[0].seVector.GetSize() ;
		struct _pair32 *chain = new struct _pair32[seCnt] ;
		bool *covered = new bool[seCnt] ;
		bool *usedConstraints = new bool[constraints.constraints.size() ] ;
		/*lowCovSubexon.Init( seCnt ) ;
		double *avgDepth = new double[seCnt ] ;

		memset( avgDepth, 0, sizeof( double ) * seCnt ) ;
		int size = constraints.constraints.size() ;
		for ( i = 0 ; i < size ; ++i )
		{
			std::vector<int> subexonIdx ;
			constraints.constraints[i].GetOnesIndices( subexonIdx ) ;
			int seIdxCnt = subexonidx.size() ;
			for ( j = 0 ; j < seIdxCnt ; ++j )
				avgDepth[ subexonidx[j] ] += constraints.constraints[i].support ;
		}
		for ( i = 0 ; i < seCnt ; ++i )
		{
			if ( avgDepth[i] * alignments.readLen / (double)( subexons[i].end - subexons[i].start + 1 ) < 1 )
		}*/
		

		for ( i = 0 ; i < atcnt ; ++i )
		{
			double value = inf ;
			int tag = -1 ;

			std::vector<int> subexonIdx ;
			alltranscripts[i].seVector.GetOnesIndices( subexonIdx ) ;
			int seIdxCnt = subexonIdx.size() ;
			for ( j = 0 ; j < seIdxCnt - 1 ; ++j )
			{
				chain[j].a = subexonIdx[j] ;
				chain[j].b = subexonIdx[j + 1] ;
				covered[j] = false ;
			}
			memset( usedConstraints, false, sizeof( bool ) * constraints.constraints.size() ) ;
			int compatibleCnt = 0 ;
			for ( j = 0 ; j < tcCnt ; ++j )
			{
				if ( btable[i].Test(j) && tc[j].abundance > 0 )
				{	
					++compatibleCnt ;
					if ( tc[j].abundance < value && 
						!( tc[j].i == tc[j].j 
							&& ( constraints.constraints[ tc[j].i ].first + constraints.constraints[ tc[j].i ].last == 2 * alltranscripts[i].first 
								||  constraints.constraints[ tc[j].i ].first + constraints.constraints[ tc[j].i ].last == 2 * alltranscripts[i].last ) ) )
					{	
						value = tc[j].abundance ;
						tag = j ;
					}
					avgTranscriptAbundance[i] += tc[j].abundance ;

					if ( !usedConstraints[ tc[j].i ] )
					{
						struct _constraint &c = constraints.constraints[ tc[j].i ] ;
						for ( k = 0 ; k < seIdxCnt - 1 ; ++k )
						{
							// Note that since the constraint is already compatible with the txpt,
							//   chain[k].a/b must be also adjacent in this constraint.
							if ( c.vector.Test( chain[k].a ) && c.vector.Test( chain[k].b ) ) 
								covered[k] = true ;
						}
						usedConstraints[ tc[j].i ] = true ;
					}

					if ( !usedConstraints[ tc[j].j ] )
					{
						struct _constraint &c = constraints.constraints[ tc[j].j ] ;
						for ( k = 0 ; k < seIdxCnt - 1 ; ++k )
						{
							if ( c.vector.Test( chain[k].a ) && c.vector.Test( chain[k].b ) )
								covered[k] = true ;
						}
						usedConstraints[ tc[j].j ] = true ;
					}
				}
			}
			// Every two-subexon chain should be covered by some reads if a transcript is expressed highly enough
			int cnt = 0 ;
			for ( j = 0 ; j < seIdxCnt - 1 ; ++j )	
				if ( covered[j] == false )
				{
					value = 0 ;
				}
				else
					++cnt ;
			if ( seIdxCnt > 1 )
				coveredPortion[i] = (double)cnt / (double)( seIdxCnt - 1 ) ;
			else
				coveredPortion[i] = 1 ;
			if ( coveredPortion[i] == 0 )
				coveredPortion[i] = (double)0.5 / ( seIdxCnt ) ;

			if ( tag == -1 ) 
				value = 0 ;
			if ( value > maxAbundance )
				maxAbundance = value ;
			transcriptAbundance[i] = value ;
			if ( tag != -1 )
				avgTranscriptAbundance[i] /= compatibleCnt ;

			//printf( "abundance %d: %lf %lf ", i, value, avgTranscriptAbundance[i] ) ;
			//alltranscripts[i].seVector.Print() ;
		}
		if ( maxAbundance == 0 )
		{
			for ( i = 0 ; i < atcnt ; ++i )
			{
				transcriptAbundance[i] = coveredPortion[i] ;
			}
			maxAbundance = 1 ;
		}
		//printf( "%s: %lf\n", __func__, maxAbundance ) ;
		delete[] usedConstraints ;
		delete[] covered ;
		delete[] chain ;
	}
	else 
	{
		for ( i = 0 ; i < atcnt ; ++i )
		{
			transcriptAbundance[i] = alltranscripts[i].abundance ;
			if ( transcriptAbundance[i] > maxAbundance )
				maxAbundance = transcriptAbundance[i] ;
			coveredPortion[i] = 1 ;
		}
		if ( maxAbundance == 0 )
			maxAbundance = 1 ;
	}
	
	int iterCnt = -1 ;
	while ( 1 )
	{
		double max = -1 ;
		int maxtag = -1 ;
		double maxcnt = -1 ;
		++iterCnt ;

		// Find the optimal candidate.
		for ( i = 0 ; i < atcnt ; ++i )
		{
			double value = inf ;
			double cnt = 0 ;
			
			for ( j = 0 ; j < tcCnt ; ++j )
			{
				if ( tc[j].abundance > 0 && btable[i].Test( j ) )
				{
					cnt += tc[j].effectiveCount ;
				}
			}

			value = transcriptAbundance[i] ;
			if ( cnt == 0 ) // This transcript does not satisfy any undepleted constraints.
				continue ;
			cnt *= coveredPortion[i] ;
			double score = ComputeScore( cnt, 1.0, value, maxAbundance, alltranscripts[i].correlationScore ) ;
			if ( cnt > maxcnt )
				maxcnt = cnt ;

			if ( score > max )
			{
				max = score ;
				maxtag = i ;
			}
			else if ( score == max )
			{
				if ( avgTranscriptAbundance[maxtag] < avgTranscriptAbundance[i] )	
				{
					max = score ;
					maxtag = i ;
				}
			}
			//printf( "score: %d %lf -> %lf\n", i, cnt, score ) ;
		}
		if ( maxcnt == 0 || maxtag == -1 )
			break ;

		// Update the abundance for each transcript.	
		double update = inf ;
		int updateTag = 0 ;
		for ( j = 0 ; j < tcCnt ; ++j )
		{
			if ( btable[ maxtag ].Test( j ) && tc[j].abundance > 0 && 
					tc[j].abundance <= update )
			{
				update = tc[j].abundance ;	
				updateTag = j ;
			}
		}

		struct _transcript nt ;
		nt.seVector.Duplicate( alltranscripts[ maxtag ].seVector ) ;
		nt.first = alltranscripts[maxtag].first ;
		nt.last = alltranscripts[maxtag].last ;
		nt.abundance = 0 ; 
		nt.partial = false ;
		for ( j = 0 ; j < tcCnt ; ++j )
		{
			if ( btable[maxtag].Test( j ) && tc[j].abundance > 0 )
			{
				tc[j].abundance -= 1 * update ;
				double factor = tc[j].effectiveCount ;
				nt.abundance += ( tc[j].support * update / tc[j].normAbund * factor ) ;
			}
			if ( tc[j].abundance < 0 )
			{
				tc[j].abundance = 0 ;
			}
		}
		//printf( "maxtag=%d %lf %lf\n", maxtag, update, nt.abundance ) ;

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

	delete[] transcriptSeCnt ;
	delete[] transcriptAbundance ;
	delete[] coveredPortion ;
}

int TranscriptDecider::RefineTranscripts( struct _subexon *subexons, std::vector<struct _transcript> &transcripts, Constraints &constraints ) 
{
	int i, j ;
	int tcnt = transcripts.size() ;
	if ( tcnt == 0 )
		return 0 ;
	int tcCnt = constraints.matePairs.size() ;

	std::vector<struct _matePairConstraint> &tc = constraints.matePairs ;
	std::vector<struct _constraint> &scc = constraints.constraints ; //single-end constraints.constraints
	
	// Remove transcripts whose FPKM are too small.
	double *geneMaxFPKM = new double[usedGeneId - baseGeneId ] ;
	memset( geneMaxFPKM, 0, sizeof( double ) * ( usedGeneId - baseGeneId ) ) ;
	double *geneMaxCov = new double[ usedGeneId - baseGeneId ] ;
	memset( geneMaxCov, 0, sizeof( double ) * ( usedGeneId - baseGeneId ) ) ;
	int *txptGid = new int[tcnt] ;
	for ( i = 0 ; i < tcnt ; ++i )
	{
		int gid = GetTranscriptGeneId( transcripts[i], subexons ) ;
		int len = GetTranscriptLengthFromAbundanceAndFPKM( transcripts[i].abundance, transcripts[i].FPKM ) ;
		//printf( "%lf %lf %d\n", transcripts[i].abundance, transcripts[i].FPKM, len ) ;
		if ( transcripts[i].FPKM > geneMaxFPKM[gid - baseGeneId ] )
			geneMaxFPKM[ gid - baseGeneId ] = transcripts[i].FPKM ;
		if ( transcripts[i].abundance * alignments.readLen / len > geneMaxCov[gid - baseGeneId ] )
			geneMaxCov[gid - baseGeneId] = ( transcripts[i].abundance * alignments.readLen ) / len ;
		txptGid[i] = gid ;
	}
	BitTable bufferTable ;
	if ( tcnt > 0 )
	{
		bufferTable.Duplicate( transcripts[0].seVector ) ;
	}

	int sccCnt = scc.size() ;
	for ( i = 0 ; i < tcnt ; ++i )
	{
		//printf( "%d: %lf %lf\n", txptGid[i], transcripts[i].abundance, geneMaxFPKM[ txptGid[i] - baseGeneId ] ) ;
		if ( transcripts[i].FPKM < FPKMFraction * geneMaxFPKM[ txptGid[i] - baseGeneId ] )
			transcripts[i].abundance = -1 ;
		if ( transcripts[i].abundance != -1 )
		{
			int len = GetTranscriptLengthFromAbundanceAndFPKM( transcripts[i].abundance, transcripts[i].FPKM ) ;
			double cov = ( transcripts[i].abundance * alignments.readLen ) / len ;
			//printf( "%d: %d %d %lf %lf\n", i, len, transcripts[i].seVector.Count(), cov, geneMaxCov[ txptGid[i] - baseGeneId ]  ) ;
			if ( ( tcnt > 1 || len <= 1000 || transcripts[i].seVector.Count() <= 3 ) && cov < txptMinReadDepth )
			{
				if ( usedGeneId == baseGeneId + 1 && transcripts[i].seVector.Count() > 3 
					&& len > 1000 && geneMaxCov[ txptGid[i] - baseGeneId ] == cov )
					continue ;
				// Test whether this transcript is fully covered. If so ,we can filter it.
				bufferTable.Reset() ;
				for ( j = 0 ; j < sccCnt ; ++j )	
				{
					if ( !IsConstraintInTranscript( transcripts[i], scc[j] ) ) 
						continue ;
					bufferTable.Or( scc[j].vector ) ;
				}
				if ( bufferTable.IsEqual( transcripts[i].seVector ) )
					transcripts[i].abundance = -1 ;
				/*else
				{
					transcripts[i].seVector.Print() ;
					bufferTable.Print() ;
					OutputTranscript( stderr, subexons, transcripts[i] ) ;
				}*/
			}
		}
	}
	
	tcnt = RemoveNegativeAbundTranscripts( transcripts )  ;
	delete []geneMaxCov ;
	bufferTable.Release() ;
	delete []geneMaxFPKM ;
	delete []txptGid ;

	/*==================================================================
	Remove transcripts that seems duplicated
	====================================================================*/
	for ( i = 0 ; i < tcnt ; ++i )
	{
		int support = 0 ;
		int uniqSupport = 0 ;

		for ( j = 0 ; j < tcCnt ; ++j )
		{
			if ( !IsConstraintInTranscript( transcripts[i], scc[ tc[j].i ] ) || !IsConstraintInTranscript( transcripts[i], scc[ tc[j].j ] ) )
				continue ;
			//support += scc[ tc[j].i ].support + scc[ tc[j].j ].support ;
			//uniqSupport += scc[ tc[j].i ].uniqSupport + scc[ tc[j].j ].uniqSupport ; 
			support += tc[j].support ;
			uniqSupport += tc[j].uniqSupport ;
		}

		if ( (double)uniqSupport < 0.05 * support )
			transcripts[i].abundance = -1 ;
	}
	tcnt = RemoveNegativeAbundTranscripts( transcripts )  ;
	
	/*==================================================================
	Remove transcripts that is too short
	====================================================================*/
	for ( i = 0 ; i < tcnt ; ++i )
	{
		int len = GetTranscriptLengthFromAbundanceAndFPKM( transcripts[i].abundance, transcripts[i].FPKM ) ;
		if ( len < 200 )
		{
			transcripts[i].abundance = -1 ;
		}
	}
	tcnt = RemoveNegativeAbundTranscripts( transcripts )  ;
	return transcripts.size() ;
}


int TranscriptDecider::Solve( struct _subexon *subexons, int seCnt, std::vector<Constraints> &constraints, SubexonCorrelation &subexonCorrelation )
{
	int i, j, k ;
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
		else if ( subexons[i].leftClassifier < canBeSoftBoundaryThreshold && subexons[i].leftClassifier != -1 
			&& subexons[i].leftStrand != 0 ) // The case of overhang.
		{
			// We then look into whether there is a left-side end already showed up before this subexon in this region of subexons.
			bool flag = true ;
			for ( j = i - 1 ; j >= 0 ; --j )
			{
				if ( subexons[j].end + 1 != subexons[j + 1].start )	
					break ;
				if ( subexons[i].canBeStart == true )
				{
					flag = false ;
					break ;
				}	
			}
			subexons[i].canBeStart = flag ;
		}

		if ( subexons[i].nextCnt == 0 )
			subexons[i].canBeEnd = true ;
		else if ( subexons[i].rightClassifier < canBeSoftBoundaryThreshold && subexons[i].rightClassifier != -1 
			&& subexons[i].rightStrand != 0 )
		{
			subexons[i].canBeEnd = true ;
		}
		// Remove other soft end already showed up in this region of subexons.
		if ( subexons[i].canBeEnd == true )
		{
			for ( j = i - 1 ; j >= 0 ; --j )
			{
				if ( subexons[j].end + 1 != subexons[j + 1].start )
					break ;
				if ( subexons[j].canBeEnd == true )
				{
					subexons[j].canBeEnd = false ;
					break ;
				}
			}
		}
		//printf( "%d: %d %lf\n", subexons[i].canBeStart, subexons[i].prevCnt, subexons[i].leftClassifier ) ;
	}

	// Go through the cases of mixture region to set canBeStart/End.
	// e.g: +[...]+_____+[....]-...]+____+[..)_____-[...]-
	//                   ^ then we need to force a start point here.
	// Do we need to associate a strand information with canBeStart, canBeEnd?
	for ( i = 0 ; i < seCnt ; )
	{
		// [i, j) is a region.
		for ( j = i + 1 ; j < seCnt ; ++j )
			if ( subexons[j].start > subexons[j - 1].end + 1 )
				break ;
		if ( subexons[i].canBeStart == false ) // then subexons[i] must has a hard left boundary.
		{
			int leftStrandCnt[2] = {0, 0} ;
			for ( k = i ; k < j ; ++k )
			{
				if ( !SubexonGraph::IsSameStrand( subexons[k].rightStrand, subexons[i].leftStrand ) )
					break ;
				if ( subexons[k].leftStrand != 0 )
					++leftStrandCnt[ ( subexons[k].leftStrand + 1 ) / 2 ] ;
			}
			if ( k < j && leftStrandCnt[ ( subexons[k].rightStrand + 1 ) / 2 ] == 0 )
				subexons[i].canBeStart = true ;
		}

		if ( subexons[j - 1].canBeEnd == false )
		{
			int rightStrandCnt[2] = {0, 0} ;
			for ( k = j - 1 ; k >= i ; --k )	
			{
				if ( !SubexonGraph::IsSameStrand( subexons[k].leftStrand, subexons[j - 1].rightStrand ) )
					break ;
				if ( subexons[k].rightStrand != 0 )
					++rightStrandCnt[ ( subexons[k].rightStrand + 1 ) / 2 ] ;
			}
			if ( k >= i && rightStrandCnt[ ( subexons[k].leftStrand + 1 ) / 2 ] == 0 )
				subexons[j - 1].canBeEnd = true ;
		}

		//if ( subexons[i].start == 6870264)
		//	printf( "hi %d %d\n",i , subexons[i].canBeStart ) ;
		i = j ;
	}
	/*for ( i = 0 ; i < seCnt ; ++i )
	{
		printf( "%d %d: %d %d\n", subexons[i].start, subexons[i].end, subexons[i].canBeStart, subexons[i].canBeEnd ) ;
	}*/

	// Find the gene ids.
	usedGeneId = baseGeneId = -1 ;
	defaultGeneId[0] = defaultGeneId[1] = -1 ;
	for ( i = 0 ; i < seCnt ; ++i )
	{
		if ( subexons[i].geneId < 0 )
			continue ;
		if ( baseGeneId == -1 || subexons[i].geneId < baseGeneId )	
			baseGeneId = subexons[i].geneId ;
		if ( subexons[i].geneId > usedGeneId )
			usedGeneId = subexons[i].geneId ;
		if ( ( subexons[i].rightStrand == -1 || subexons[i].leftStrand == -1 ) &&
			( defaultGeneId[0] == -1 || subexons[i].geneId < defaultGeneId[0] ) )
			defaultGeneId[0] = subexons[i].geneId ;
		if ( ( subexons[i].rightStrand == 1 || subexons[i].leftStrand == 1 ) &&
			( defaultGeneId[1] == -1 || subexons[i].geneId < defaultGeneId[1] ) )
			defaultGeneId[1] = subexons[i].geneId ;
			defaultGeneId[0] = subexons[i].geneId ;
	}
	++usedGeneId ;
	//printf( "%d %d\n", baseGeneId, usedGeneId ) ;
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
		{
			double msize = constraints[i].matePairs.size() ;
			double csize = constraints[i].constraints.size() ;
			if ( cnt > ( csize / msize ) * ( csize / msize ) * seCnt 
				&& cnt > USE_DP / ( msize * msize ) && cnt > 50 )
			{
				useDP = true ;
				break ;
			}
		}
	}

	int atCnt = cnt ;
	printf( "atCnt=%d %d %d %d\n", atCnt, useDP, (int)constraints[0].constraints.size(), (int)constraints[0].matePairs.size() ) ;
	std::vector<struct _transcript> alltranscripts ;
	//useDP = false ;	
	if ( !useDP )
	{
		alltranscripts.resize( atCnt ) ;
		for ( i = 0 ; i < atCnt ; ++i )
		{
			alltranscripts[i].seVector.Init( seCnt ) ; 
			alltranscripts[i].correlationScore = 1 ;
		}

		atCnt = 0 ;
		for ( i = 0 ; i < seCnt ; ++i )
		{
			if ( subexons[i].canBeStart )
				EnumerateTranscript( i, 0, f, 0, subexons, subexonCorrelation, 1, alltranscripts, atCnt ) ;
		}
		alltranscripts.resize( atCnt ) ;
		//printf( "transcript cnt: %d\n", atCnt ) ;
		//printf( "%d %d\n", alltranscripts[0].seVector.Test( 1 ), constraints[0].matePairs.size() ) ;
	}
	else // Use dynamic programming to pick a set of candidate transcript.
	{
		std::vector<struct _transcript> sampleTranscripts ;

		// pre allocate the memory.
		struct _dpAttribute attr ;
		attr.f1 = new struct _dp[seCnt] ;
		attr.f2 = new struct _dp*[seCnt] ;
		for ( i = 0 ; i < seCnt ; ++i )
			attr.f2[i] = new struct _dp[seCnt] ;
		attr.hash = new struct _dp[HASH_MAX] ;
		for ( i = 0 ; i < seCnt ; ++i )
		{
			attr.f1[i].seVector.Nullify() ;
			attr.f1[i].seVector.Init( seCnt ) ;
			for ( j = i ; j < seCnt ; ++j )
			{
				attr.f2[i][j].seVector.Nullify() ;
				attr.f2[i][j].seVector.Init( seCnt ) ;
			}
		}
		for ( i = 0 ; i < HASH_MAX ; ++i )
		{
			attr.hash[i].seVector.Nullify() ;
			attr.hash[i].seVector.Init( seCnt ) ;
		}

		// select candidate transcripts from each sample.
		for ( i = 0 ; i < sampleCnt ; ++i )
		{
			sampleTranscripts.clear() ;
			PickTranscriptsByDP( subexons, seCnt, constraints[i], subexonCorrelation, attr, sampleTranscripts ) ;		
			int size = sampleTranscripts.size() ;
			for ( j = 0 ; j < size ; ++j )
				alltranscripts.push_back( sampleTranscripts[j] ) ;
		}

		// release the memory.
		for ( i = 0 ; i < seCnt ; ++i )	
		{
			attr.f1[i].seVector.Release() ;
			for ( j = i ; j < seCnt ; ++j )
				attr.f2[i][j].seVector.Release() ;
		}
		for ( i = 0 ; i < HASH_MAX ; ++i )
			attr.hash[i].seVector.Release() ;

		delete[] attr.f1 ;
		for ( i = 0 ; i < seCnt ; ++i )
			delete[] attr.f2[i] ;
		delete[] attr.f2 ;
		delete[] attr.hash ;

		// we can further pick a smaller subsets of transcripts here if the number is still to big. 
		CoalesceSameTranscripts( alltranscripts ) ;
	}

	transcriptId = new int[usedGeneId - baseGeneId] ;
	for ( i = 0 ; i < sampleCnt ; ++i )
	{
		std::vector<struct _transcript> predTranscripts ;
		int size = alltranscripts.size() ;
		for ( j = 0 ; j < size ; ++j )
			alltranscripts[j].abundance = -1 ;
		PickTranscripts( alltranscripts, constraints[i], subexonCorrelation, predTranscripts ) ;

		size = predTranscripts.size() ;
		InitTranscriptId() ;
		for ( j = 0 ; j < size ; ++j )
		{
			ConvertTranscriptAbundanceToFPKM( subexons, predTranscripts[j] ) ;
		}
		size = RefineTranscripts( subexons, predTranscripts, constraints[i] ) ;
		for ( j = 0 ; j < size ; ++j )
		{
			OutputTranscript( i, subexons, predTranscripts[j] ) ;
		}
		for ( j = 0 ; j < size ; ++j )
			predTranscripts[j].seVector.Release() ;
	}
	delete []transcriptId ;

	atCnt = alltranscripts.size() ;
	for ( i = 0 ; i < atCnt ; ++i )
		alltranscripts[i].seVector.Release() ;
	compatibleTestVectorT.Release() ;
	compatibleTestVectorC.Release() ;
	delete[] f ;
	return 0 ;	
}

void *TranscriptDeciderSolve_Wrapper( void *a ) 
{
	int i ;

	struct _transcriptDeciderThreadArg &arg = *( (struct _transcriptDeciderThreadArg *)a ) ;
	TranscriptDecider transcriptDecider( arg.FPKMFraction, arg.classifierThreshold, arg.txptMinReadDepth, arg.sampleCnt, *( arg.alignments ) ) ;
	transcriptDecider.SetNumThreads( arg.numThreads + 1 ) ;
	transcriptDecider.SetMultiThreadOutputHandler( arg.outputHandler ) ;
	transcriptDecider.Solve( arg.subexons, arg.seCnt, arg.constraints, arg.subexonCorrelation ) ;
	
	// Release memory
	for ( i = 0 ; i < arg.seCnt ; ++i )
	{
		delete[] arg.subexons[i].prev ;
		delete[] arg.subexons[i].next ;
	}
	delete[] arg.subexons ;

	// Put the work id back to the free threads queue.
	pthread_mutex_lock( arg.ftLock ) ;
	arg.freeThreads[ *( arg.ftCnt ) ] = arg.tid ;
	++*( arg.ftCnt ) ;
	if ( *( arg.ftCnt ) == 1 )
		pthread_cond_signal( arg.fullWorkCond ) ;
	pthread_mutex_unlock( arg.ftLock) ;

	pthread_exit( NULL ) ;
}

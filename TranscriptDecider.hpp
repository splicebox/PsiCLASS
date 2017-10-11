#ifndef _MOURISL_CLASSES_TRANSCRIPTDECIDER_HEADER
#define _MOURISL_CLASSES_TRANSCRIPTDECIDER_HEADER

#include "alignments.hpp"
#include "SubexonGraph.hpp"
#include "SubexonCorrelation.hpp"
#include "BitTable.hpp"
#include "Constraints.hpp"

struct _transcript
{
	BitTable seVector ;
	double abundance ;
	double correlationScore ; 
	double FPKM ;

	int first, last ; // indicate the index of the first and last subexons.
	bool partial ; // wehther this is a partial transcript.
} ;

struct _dp
{
	BitTable seVector ; 
	int first, last ;
	// The "cnt" is for the hash structure.
	// the first cnt set bits represent the subexons that are the key of the hash
	// the remaining set bits are the optimal subtranscript follow the key.
	int cnt ; 
	double cover ;

	double minAbundance ;
	int timeStamp ;
	int strand ;
} ;


struct _dpAttribute
{
	struct _dp *f1, **f2 ;
	struct _dp *hash ;

	struct _transcript bufferTxpt ;

	bool forAbundance ;

	struct _subexon *subexons ;
	int seCnt ;

	double minAbundance ;
	int timeStamp ;
} ;

class TranscriptDecider
{
private:
	int sampleCnt ;
	double FPKMFraction ;

	Constraints *constraints ;
	//struct _subexon *subexons ;
	//int seCnt ;

	int *geneId ; // assign the gene id to each subexon in this region.
	int usedGeneId ;
	int baseGeneId, defaultGeneId[2] ;

	int *transcriptId ; // the next transcript id for each gene id (we shift the gene id to 0 in this array.)
	Alignments &alignments ; // for obtain the chromosome names.

	std::vector<FILE *> outputFPs ;

	BitTable compatibleTestVectorT, compatibleTestVectorC ;

	// Test whether subexon tag is a start subexon in a mixture region that corresponds to the start of a gene on another strand.
	bool IsStartOfMixtureStrandRegion( int tag, struct _subexon *subexons, int seCnt ) ;

	// The functions to pick transcripts through dynamic programming
	void SearchSubTranscript( int tag, int strand, int parents[], int pcnt, struct _dp &pdp, int visit[], int vcnt, int extends[], int extendCnt, std::vector<struct _constraint> &tc, int tcStartInd, struct _dpAttribute &attr ) ;
	struct _dp SolveSubTranscript( int visit[], int vcnt, int strand, std::vector<struct _constraint> &tc, int tcStartInd, struct _dpAttribute &attr ) ;
	void PickTranscriptsByDP( struct _subexon *subexons, int seCnt, Constraints &constraints, SubexonCorrelation &correlation, struct _dpAttribute &attr, std::vector<struct _transcript> &allTranscripts ) ;

	void SetDpContent( struct _dp &a, struct _dp &b, const struct _dpAttribute &attr )
	{
		a.seVector.Assign( b.seVector ) ;
		a.first = b.first ;
		a.last = b.last ;
		a.cnt = b.cnt ;
		a.cover = b.cover ;
		
		a.strand = b.strand ;
		a.minAbundance = attr.minAbundance ;
		a.timeStamp = attr.timeStamp ;
	} 

	void ResetDpContent( struct _dp &d )
	{
		d.seVector.Reset() ;
		d.first = -1 ;
		d.last = -1 ;
		d.cnt = -1 ;
		d.cover = -1 ;
		d.minAbundance = -1 ;
		d.timeStamp = -1 ;
	}

	double canBeSoftBoundaryThreshold ;

	// Test whether a constraints is compatible with the transcript.
	// Return 0 - uncompatible or does not overlap at all. 1 - fully compatible. 2 - Head of the constraints compatible with the tail of the transcript
	int IsConstraintInTranscript( struct _transcript transcript, struct _constraint &c ) ;
	int IsConstraintInTranscriptDebug( struct _transcript transcript, struct _constraint &c ) ;
	
	// Count how many transcripts are possible starting from subexons[tag].
	int SubTranscriptCount( int tag, struct _subexon *subexons, int f[] ) ;

	// The methods when there is no need for DP
	void EnumerateTranscript( int tag, int strand, int visit[], int vcnt, struct _subexon *subexons, SubexonCorrelation &correlation, double correlationScore, std::vector<struct _transcript> &alltranscripts, int &atcnt ) ;
	// For the simpler case, we can pick sample by sample.
	void PickTranscripts( std::vector<struct _transcript> &alltranscripts, Constraints &constraints, SubexonCorrelation &seCorrelation, std::vector<struct _transcript> &transcripts ) ; 

	static bool CompSortTranscripts( const struct _transcript &a, const struct _transcript &b )
	{
		if ( a.first < b.first )
			return true ;
		else if ( a.first > b.first )
			return false ;

		int diffPos = a.seVector.GetFirstDifference( b.seVector ) ;
		if ( diffPos == -1 )
			return false ;
		if ( a.seVector.Test( diffPos ) )
			return true ;
		else
			return false ;
	} 

	static int CompExtendsPairs( const void *p1, const void *p2 )
	{
		return ((struct _pair32 *)p1)->b - ((struct _pair32 *)p2)->b ;
	}

	double ComputeScore( double cnt, double seCnt, double a, double A, double correlation )
	{
		return ( cnt / seCnt ) * ( 1 + pow( a / A, 0.25 ) ) + correlation ;
		//return ( cnt ) * ( exp( 1 + a / A ) ) + correlation ; 
	}


	void ConvertTranscriptAbundanceToFPKM( struct _subexon *subexons, struct _transcript &t )
	{
		int txptLen = 0 ;
		int i, size ;

		std::vector<int> subexonInd ;
		t.seVector.GetOnesIndices( subexonInd ) ;
		size = subexonInd.size() ;
		for ( i = 0 ; i < size ; ++i )
			txptLen += ( subexons[ subexonInd[i] ].end - subexons[ subexonInd[i] ].start + 1 ) ;
		t.FPKM = t.abundance / ( ( alignments.totalReadCnt / 1000000.0 ) * ( txptLen / 1000.0 ) ) ;
	}

	void CoalesceSameTranscripts( std::vector<struct _transcript> &t ) ;

	// The function to assign gene ids to subexons.
	void SetGeneId( int tag, int strand, struct _subexon *subexons, int seCnt, int id ) ;

	// Initialize the structure to store transcript id 
	void InitTranscriptId() ; 
	
	int GetTranscriptGeneId( std::vector<int> &subexonInd, struct _subexon *subexons ) ;
	int GetTranscriptGeneId( struct _transcript &t, struct _subexon *subexons ) ;
	int RemoveNegativeAbundTranscripts( std::vector<struct _transcript> &transcripts )
	{
		int i, j ;
		int tcnt = transcripts.size() ;
		j = 0 ;
		for ( i = 0 ; i < tcnt ; ++i )
		{
			if ( transcripts[i].abundance == -1 )
			{
				transcripts[i].seVector.Release() ; // Don't forget release the memory.
				continue ;
			}
			transcripts[j] = transcripts[i] ;
			++j ;
		}
		transcripts.resize( j ) ;
		return j ;
	}

	int RefineTranscripts( struct _subexon *subexons, std::vector<struct _transcript> &transcripts, Constraints &constraints ) ;

	void OutputTranscript( FILE *fp, struct _subexon *subexons, struct _transcript &transcript ) ;
public:
	TranscriptDecider( double f, double c, int sampleCnt, Alignments &a ): alignments( a )  
	{
		FPKMFraction = f ;
		canBeSoftBoundaryThreshold = c ;
		usedGeneId = 0 ;
		defaultGeneId[0] = -1 ;
		defaultGeneId[1] = -1 ;
		this->sampleCnt = sampleCnt ;
	}
	~TranscriptDecider() 
	{
		int i ;
		int size = outputFPs.size() ;
		for ( i = 0 ; i < size ; ++i )
		{
			fclose( outputFPs[i] ) ;
		}
	}
	
	// @return: the number of assembled transcript 
	int Solve( struct _subexon *subexons, int seCnt, std::vector<Constraints> &constraints, SubexonCorrelation &subexonCorrelation ) ;

	void SetOutputFPs( char *outputPrefix )
	{
		int i ;
		char buffer[1024] ;
		for ( i = 0 ; i < sampleCnt ; ++i )
		{
			if ( outputPrefix[0] )
				sprintf( buffer, "%s_sample_%d.gtf", outputPrefix, i ) ;
			else
				sprintf( buffer, "sample_%d.gtf", i ) ;
			FILE *fp = fopen( buffer, "w" ) ;
			outputFPs.push_back( fp ) ;
		}
	}
} ;

#endif

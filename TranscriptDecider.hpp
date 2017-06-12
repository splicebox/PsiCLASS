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

	int first, last ; // indicate the index of the first and last subexons.
	bool partial ; // wehther this is a partial transcript.
} ;


class TranscriptDecider
{
private:
	int sampleCnt ;
	
	Constraints *constraints ;
	struct _subexon *subexons ;
	int seCnt ;

	BitTable compatibleTestVectorT, compatibleTestVectorC ;

	void SolveByEnumeration() ;
	void PickTranscriptByDP() ;

	double canBeSoftBoundaryThreshold ;

	// Test whether a constraints is compatible with the transcript.
	// Return 0 - uncompatible or does not overlap at all. 1 - fully compatible. 2 - Head of the constraints compatible with the tail of the transcript
	int IsConstraintInTranscript( struct _transcript transcript, struct _constraint &c ) ;
	
	// Count how many transcripts are possible starting from subexons[tag].
	int SubTranscriptCount( int tag, struct _subexon *subexons, int f[] ) ;

	// The methods when there is no need for DP
	void EnumerateTranscript( int tag, int visit[], int vcnt, struct _subexon *subexons, SubexonCorrelation &correlation, double correlationScore, struct _transcript *alltranscripts, int &atcnt ) ;
	// For the simpler case, we can pick sample by sample.
	void PickTranscripts( struct _transcript *alltranscripts, const int &atcnt, Constraints &constraints, SubexonCorrelation &seCorrelation, std::vector<struct _transcript> &transcripts ) ; 

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

	double ComputeScore( double cnt, double a, double A, double correlation )
	{
		return cnt * ( 1 + pow( a / A, 0.25 ) ) + correlation ;
	}

	void CoalesceSameTranscripts( std::vector<struct _transcript> &t ) ;

public:
	TranscriptDecider( int sampleCnt ) 
	{
		canBeSoftBoundaryThreshold = 0.05 ;
		this->sampleCnt = sampleCnt ;
	}
	~TranscriptDecider() {}
	
	// @return: the number of assembled transcript 
	int Solve( struct _subexon *subexons, int seCnt, Constraints constraints[], SubexonCorrelation &subexonCorrelation ) ;
} ;

#endif

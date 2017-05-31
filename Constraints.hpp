#ifndef _MOURISL_CLASSES_CONSTRAINTS_HEADER
#define _MOURISL_CLASSES_CONSTRAINTS_HEADER

#include <vector> 
#include <algorithm>

#include "BitTable.hpp"
#include "alignments.hpp"
#include "SubexonGraph.hpp"

struct _constraint
{
	BitTable vector ; // subexon vector
	double weight ;
	double abund ;
	double normAbund ;
	int support ;
} ;

struct _mateIdx
{
	int i, j ;
	int support ;
	double abund ;
} ;

class Constraints
{
private:
	std::vector<struct _constraint> constraints ;
	int prevStart, prevEnd ;	
	std::vector<struct _mateIdx> mates ; 

	Alignments &alignments ;
	
	//@return: whether this alignment is compatible with the subexons or not.
	bool ConvertAlignmentToBitTable( struct _pair *segments, int segCnt, struct _subexon *subexons, int seCnt, int seStart, struct _constraint &ct ) ;

	// Sort to increasing order. Since the first subexon occupies the least important digit.
	static bool CompSortConstraints( const struct _constraint &a, const struct _constraint &b )
	{
		//int k 
		int diffPos = a.vector.GetFirstDifference( b.vector ) ;
		if ( diffPos == -1 )
			return false ;

		if ( a.vector.Test( diffPos ))
			return true ;
		else
			return false ;
	}
	void CoalesceSameConstraints() ;
public:
	Constraints( Alignments &a ): alignments( a ) 
	{
	}
	~Constraints() 
	{
	}

	void Clear() 
	{
		//TODO: do I need to release the memory from BitTable?
		constraints.clear() ;
	}

	int BuildConstraints( struct _subexon *subexons, int scnt, int start, int end ) ;
} ;

#endif

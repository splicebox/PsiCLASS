#include <stdio.h>
#include <getopt.h>
#include <vector>

#include "alignments.hpp"
#include "SubexonGraph.hpp"
#include "SubexonCorrelation.hpp"
#include "Constraints.hpp"

char usage[] = "./classes [OPTIONS]:\n"
	"Required:\n"
	"\t-s STRING: path to the subexon file.\n"
	"\t-b STRING: path to the BAM file.\n"
	"Optional:\n"
	"\t--sl STRING: path to the list of single-sample subexon files.\n"
	"\t-c FLOAT: only use the subexons with classifier score <= than the given number. (default: 0.05)\n" ;

static const char *short_options = "s:b:h" ;
static struct option long_options[] =
	{
		{ "sl", required_argument, 0, 10000 },
		{ (char *)0, 0, 0, 0} 
	} ;


double classifierThreshold ;

int main( int argc, char *argv[] )
{
	int i, j, k ;

	if ( argc <= 1 )
	{
		printf( "%s", usage ) ;
		return 0 ;
	}

	int c, option_index ; // For getopt
	option_index = 0 ;
	FILE *fpSubexon = NULL ;
	std::vector<Alignments> alignmentFiles ;
	SubexonCorrelation subexonCorrelation ;

	while ( 1 )
	{
		c = getopt_long( argc, argv, short_options, long_options, &option_index ) ;
		if ( c == -1 )
			break ;

		if ( c == 's' )
		{
			fpSubexon = fopen( optarg, "r" ) ;
		}
		else if ( c == 'b' )
		{
			Alignments a ;
			a.Open( optarg ) ;
			alignmentFiles.push_back( a ) ;
		}
		else if ( c == 10000 )
		{
			subexonCorrelation.Initialize( optarg ) ;
		}
		else
		{
			printf( "%s", usage ) ;
			exit( 1 ) ;
		}
	}
	if ( fpSubexon == NULL )			
	{
		printf( "Must use -s option to speicfy subexon file.\n" ) ;
		exit( 1 ) ;
	}
	if ( alignmentFiles.size() < 1 )
	{
		printf( "Must use -b option to specify BAM files.\n" ) ;
		exit( 1 ) ;
	}

	// Build the subexon graph
	SubexonGraph subexonGraph( classifierThreshold, alignmentFiles[0], fpSubexon ) ;
	subexonGraph.ComputeGeneIntervals() ;
	
	// Solve gene by gene
	int giCnt = subexonGraph.geneIntervals.size() ;
	Constraints constraints( alignmentFiles[0] ) ;
	for ( i = 0 ; i < giCnt ; ++i )
	{
		struct _geneInterval gi = subexonGraph.geneIntervals[i] ;
		printf( "%d: %d %d %d\n", i, gi.endIdx - gi.startIdx + 1, gi.start, gi.end ) ;	
		struct _subexon *intervalSubexons = new struct _subexon[ gi.endIdx - gi.startIdx + 1 ] ;
		subexonGraph.ExtractSubexons( gi.startIdx, gi.endIdx, intervalSubexons ) ;
		
		subexonCorrelation.ComputeCorrelation( intervalSubexons, gi.endIdx - gi.startIdx + 1, alignmentFiles[0] ) ;		
		constraints.BuildConstraints( intervalSubexons, gi.endIdx - gi.startIdx + 1, gi.start, gi.end ) ;	
		
		delete[] intervalSubexons ;
	}
	
	//
	return 0 ;
}

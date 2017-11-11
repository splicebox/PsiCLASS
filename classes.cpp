#include <stdio.h>
#include <getopt.h>
#include <vector>

#include "alignments.hpp"
#include "SubexonGraph.hpp"
#include "SubexonCorrelation.hpp"
#include "Constraints.hpp"
#include "TranscriptDecider.hpp"

char usage[] = "./classes [OPTIONS]:\n"
	"Required:\n"
	"\t-s STRING: path to the subexon file.\n"
	"\t-b STRING: path to the BAM file.\n"
	"\t\tor\n"
	"\t--lb STRING: path to the file of the list of BAM files.\n"
	"Optional:\n"
	"\t--ls STRING: path to the file of the list of single-sample subexon files. (default: not used)\n"
	"\t-o STRING: the prefix of the output file. (default: not used)"
	"\t-c FLOAT: only use the subexons with classifier score <= than the given number. (default: 0.05)\n" 
	"\t-f FLOAT: filter the transcript from the gene if its abundance is lower than the given number percent of the most abundant one. (default: 0.01)\n"
	"\t-d FLOAT: filter the transcript whose average read depth is less than the given number. (default: 2.5)\n"
	;

static const char *short_options = "s:b:f:o:d:h" ;
static struct option long_options[] =
	{
		{ "ls", required_argument, 0, 10000 },
		{ "lb", required_argument, 0, 10001 },
		{ (char *)0, 0, 0, 0} 
	} ;



int main( int argc, char *argv[] )
{
	int i, j ;
	int size ;

	if ( argc <= 1 )
	{
		printf( "%s", usage ) ;
		return 0 ;
	}

	int c, option_index ; // For getopt
	option_index = 0 ;
	FILE *fpSubexon = NULL ;
	double FPKMFraction = 0.01 ; 
	double classifierThreshold ;
	double txptMinReadDepth = 2.5 ;
	char outputPrefix[1024] = "" ;
	
	std::vector<Alignments> alignmentFiles ;
	SubexonCorrelation subexonCorrelation ;
	
	classifierThreshold = 0.05 ;
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
		else if ( c == 'f' )
		{
			FPKMFraction = atof( optarg ) ; 		
		}
		else if ( c == 'o' )
		{
			strcpy( outputPrefix, optarg ) ;	
		}
		else if ( c == 'd' )
		{
			txptMinReadDepth = atof( optarg ) ;
		}
		else if ( c == 10000 ) // the list of subexon files.
		{
			subexonCorrelation.Initialize( optarg ) ;
		}
		else if ( c == 10001 ) // the list of bam files.
		{
			FILE *fp = fopen( optarg, "r" ) ;
			char buffer[1024] ;
			while ( fgets( buffer, sizeof( buffer ), fp ) != NULL )
			{
				int len = strlen( buffer ) ;
				if ( buffer[len - 1] == '\n' )
				{
					buffer[len - 1] = '\0' ;
					--len ;

				}
				Alignments a ;
				a.Open( buffer ) ;
				alignmentFiles.push_back( a ) ;
			}
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
	
	size = alignmentFiles.size() ;
	for ( i = 0 ; i < size ; ++i )
	{
		alignmentFiles[i].GetGeneralInfo() ;
		alignmentFiles[i].Rewind() ;
	}

	// Build the subexon graph
	SubexonGraph subexonGraph( classifierThreshold, alignmentFiles[0], fpSubexon ) ;
	subexonGraph.ComputeGeneIntervals() ;
	
	// Solve gene by gene
	int sampleCnt = alignmentFiles.size() ;
	std::vector<Constraints> multiSampleConstraints ;
	for ( i = 0 ; i < sampleCnt ; ++i )
	{
		Constraints constraints( &alignmentFiles[i] ) ;
		multiSampleConstraints.push_back( constraints ) ;
	}
	TranscriptDecider transcriptDecider( FPKMFraction, classifierThreshold, txptMinReadDepth, sampleCnt, alignmentFiles[0] ) ;

	transcriptDecider.SetOutputFPs( outputPrefix ) ;
	
	int giCnt = subexonGraph.geneIntervals.size() ;
	for ( i = 0 ; i < giCnt ; ++i )
	{
		struct _geneInterval gi = subexonGraph.geneIntervals[i] ;
		printf( "%d: %d %d %d\n", i, gi.endIdx - gi.startIdx + 1, gi.start, gi.end ) ;	
		struct _subexon *intervalSubexons = new struct _subexon[ gi.endIdx - gi.startIdx + 1 ] ;
		subexonGraph.ExtractSubexons( gi.startIdx, gi.endIdx, intervalSubexons ) ;
		
		subexonCorrelation.ComputeCorrelation( intervalSubexons, gi.endIdx - gi.startIdx + 1, alignmentFiles[0] ) ;
		for ( j = 0 ; j < sampleCnt ; ++j )
			multiSampleConstraints[j].BuildConstraints( intervalSubexons, gi.endIdx - gi.startIdx + 1, gi.start, gi.end ) ;	
		
		transcriptDecider.Solve( intervalSubexons, gi.endIdx - gi.startIdx + 1, multiSampleConstraints, subexonCorrelation ) ;

		for ( j = 0 ; j < gi.endIdx - gi.startIdx + 1 ; ++j )
		{
			delete[] intervalSubexons[j].prev ;
			delete[] intervalSubexons[j].next ;
		}
		delete[] intervalSubexons ;
	}
	
	//
	return 0 ;
}

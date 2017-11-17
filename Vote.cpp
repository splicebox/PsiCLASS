// The program that vote to pick the trusted transcripts.
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <algorithm>
#include <string.h>
#include <string>

#include "defs.h"
#include "TranscriptDecider.hpp"

char usage[] = "./transcript-vote [OPTIONS] > output.gtf:\n"
	"Required:\n"
	"\t--lg: path to the list of GTF files.\n"
	"Optional:\n" 
	"\t-f FLOAT: the fraction of samples the transcript showed up. (default: 0.1)\n"
	"\t-n INT: the number of samples a transcript showed up. (default: 3)\n"
	;


/*struct _transcript
{
	int chrId ;
	int geneId ; 
	char strand ;
	struct _pair32 *exons ;
	int ecnt ;
	int sampleId ;
} ;*/

int GetTailNumber( char *s )
{
	int len = strlen( s ) ;
	int ret = 0 ;
	int i ;
	int factor = 1 ;
	for ( i = len - 1 ; i >= 0 && s[i] >= '0' && s[i] <= '9' ; --i, factor *= 10 )
	{
		ret += factor * ( s[i] - '0' ) ;
	}
	return ret ;
}

int CompTranscripts( const struct _outputTranscript &a, const struct _outputTranscript &b )
{
	int i ;
	if ( a.geneId != b.geneId )
		return a.geneId - b.geneId ;
	if ( a.ecnt != b.ecnt )
		return a.ecnt - b.ecnt ;

	for ( i = 0 ; i < a.ecnt ; ++i )
	{
		if ( a.exons[i].a != b.exons[i].a )					
			return a.exons[i].a - b.exons[i].a ;

		if ( a.exons[i].b != b.exons[i].b )					
			return a.exons[i].b - b.exons[i].b ;
	}
	return 0 ;
}

bool CompSortTranscripts( const struct _outputTranscript &a, const struct _outputTranscript &b )
{
	int tmp = CompTranscripts( a, b ) ;
	if ( tmp < 0 )
		return true ;
	else if ( tmp > 0 )
		return false ;
	else
		return a.sampleId < b.sampleId ;
}

int main( int argc, char *argv[] )
{
	int i, j ;
	int fraction = 0.1 ;
	int minSampleCnt = 3 ;
	std::map<std::string, int> chrNameToId ;
	std::map<int, std::string> chrIdToName ;

	FILE *fpGTFlist = NULL ;
	FILE *fp = NULL ;
	for ( i = 1 ; i < argc ; ++i )
	{
		if ( !strcmp( argv[i], "--lg" ) )
		{
			fpGTFlist = fopen( argv[i + 1], "r" ) ;
			++i ;
		}
		else if ( !strcmp( argv[i], "-f" ) )
		{
			fraction = atof( argv[i + 1] ) ;
			++i ;
		}
		else if ( !strcmp( argv[i], "-c" ) )
		{
			minSampleCnt = atoi( argv[i + 1] ) ;
			++i ;
		}
		else 
		{
			printf( "%s", usage ) ;
			exit( 1 ) ;
		}
	}
	if ( fpGTFlist == NULL )
	{
		printf( "Must use --lg option to speicfy the list of GTF files.\n" ) ;
		exit( 1 ) ;
	}
	
	std::vector<struct _outputTranscript> transcripts ;
	std::vector<struct _outputTranscript> outputTranscripts ;

	char buffer[4096] ;
	char line[10000] ;
	char chrom[50], tool[20], type[40], strand[3] ;
	//char tid[50] ;
	int start, end ;
	std::vector<struct _pair32> tmpExons ;
	int sampleCnt = 0 ;
	int gid = -1 ;
	int chrId = -1 ;
	int chrIdUsed = 0 ;
	char cStrand ;
	char prefix[50] ;

	while ( fgets( buffer, sizeof( buffer ), fpGTFlist ) != NULL )
	{
		int len = strlen( buffer ) ;	
		if ( buffer[len - 1] == '\n' )
		{
			buffer[len - 1] = '\0' ; 
			--len ;
		}
		fp = fopen( buffer, "r" ) ;
		while ( fgets( line, sizeof( line ), fp ) != NULL )
		{
			sscanf( line, "%s %s %s %d %d %s %s", chrom, tool, type, &start, &end, buffer, strand ) ;

			if ( strcmp( type, "exon" ) )
			{
				if ( tmpExons.size() > 0 )
				{
					struct _outputTranscript nt ;
					nt.sampleId = sampleCnt ; 
					nt.chrId = chrId ;
					nt.geneId = gid ;
					nt.strand = cStrand ;
					nt.ecnt = tmpExons.size() ;
					nt.exons = new struct _pair32[nt.ecnt] ;
					for ( i = 0 ; i < nt.ecnt ; ++i )
						nt.exons[i] = tmpExons[i] ;
					tmpExons.clear() ;
					transcripts.push_back( nt ) ;
				}
				continue ;
			}

			if ( chrNameToId.find( std::string( chrom ) ) == chrNameToId.end() )
			{
				chrId = chrIdUsed ;

				std::string s( chrom ) ;
				chrNameToId[s] = chrIdUsed ;
				chrIdToName[ chrIdUsed ] = s ;
				++chrIdUsed ;
			}
			else
				chrId = chrNameToId[ std::string( chrom ) ] ;

			char *p = strstr( line, "gene_id" ) ;
			for ( ; *p != ' ' ; ++p )
				;
			p += 2 ; // add extra 1 to skip \"
			sscanf( p, "%s", buffer ) ;
			//printf( "+%s %d\n", tid, strlen( tid )  ) ;
			p = buffer + strlen( buffer ) ;
			while ( *p != '\"' )
				--p ;
			*p = '\0' ;
			gid = GetTailNumber( buffer ) ;

			cStrand = strand[0] ;

			// Look for the user-defined prefix.
			int len = strlen( buffer ) ;
			j = 0 ;
			for ( i = len - 1 ; i >= 0 ; --i )
			{
				if ( buffer[i] == '.' )
				{
					++j ;
					if ( j >= 2 )
						break ;
				}
			}

			for ( j = 0 ; j < i ; ++i )
			{
				prefix[j] = buffer[j] ;
			}
			if ( i > 0 )
			{
				prefix[j] = '.' ;
				prefix[j + 1] = '\0' ;
			}
			else
				prefix[0] = '\0' ;

			struct _pair32 ne ;
			ne.a = start ;
			ne.b = end ;
			tmpExons.push_back( ne ) ;
		}
		if ( tmpExons.size() > 0 )
		{
			struct _outputTranscript nt ;
			nt.sampleId = sampleCnt ; 
			nt.chrId = chrId ;
			nt.geneId = gid ;
			nt.strand = cStrand ;
			nt.ecnt = tmpExons.size() ;
			nt.exons = new struct _pair32[nt.ecnt] ;
			for ( i = 0 ; i < nt.ecnt ; ++i )
				nt.exons[i] = tmpExons[i] ;
			tmpExons.clear() ;
			transcripts.push_back( nt ) ;
		}

		++sampleCnt ;
		fclose( fp ) ;
	}
	fclose( fpGTFlist ) ;
	std::sort( transcripts.begin(), transcripts.end(), CompSortTranscripts ) ;
	int size = transcripts.size() ;
	for ( i = 0 ; i < size ; )
	{
		for ( j = i + 1 ; j < size ; ++j )
		{
			if ( CompTranscripts( transcripts[i], transcripts[j] ) )
				break ;
		}
		// [i,j) are the same transcripts.
		if ( j - i >= fraction * sampleCnt && j - i >= minSampleCnt )
		{
			outputTranscripts.push_back( transcripts[i] ) ;
		}
		i = j ;
	}
	
	size = outputTranscripts.size() ;
	printf( "%d\n", size ) ;
	int transcriptId = 0 ;
	int prevGid = -1 ;
	for ( i = 0 ; i < size ; ++i )
	{
		struct _outputTranscript &t = outputTranscripts[i] ;
		const char *chrom = chrIdToName[t.chrId].c_str() ;
		if ( t.geneId != prevGid )
			transcriptId = 0 ;
		else
			++transcriptId ;

		fprintf( stdout, "%s\tCLASSES\ttranscript\t%d\t%d\t1000\t%c\t.\tgene_id \"%s%s.%d\"; transcript_id \"%s%s.%d.%d\";\n",
				chrom, t.exons[0].a, t.exons[t.ecnt - 1].b, t.strand,
				prefix, chrom, t.geneId,
				prefix, chrom, t.geneId, transcriptId ) ;
		for ( j = 0 ; j < t.ecnt ; ++j )
		{
			fprintf( stdout, "%s\tCLASSES\texon\t%d\t%d\t1000\t%c\t.\tgene_id \"%s%s.%d\"; "
					"transcript_id \"%s%s.%d.%d\"; exon_number \"%d\";\n",
					chrom, t.exons[j].a, t.exons[j].b, t.strand,
					prefix, chrom, t.geneId,
					prefix, chrom, t.geneId, transcriptId,
					j + 1 ) ;
		}
		prevGid = t.geneId ;
	}
	return 0 ;
}

#include <stdio.h>
#include <string.h>

#include <algorithm>
#include <vector>
#include <map>
#include <string>

#include "alignments.hpp"
#include "blocks.hpp"
#include "stats.hpp"

char usage[] = "combineSubexons [options]\n"
	       "Required options:\n"
	       "\t-s STRING: the path to the predicted subexon information. Can use multiple -s to specify multiple subexon prediction files\n" ;

struct _overhang
{
	int cnt ; //  the number of samples support this subexon.
	int length ;
	double classifier ;
} ;

struct _intronicInfo
{
	int chrId ;
	int start, end ;
	double irClassifier ;
	int irCnt ;
	struct _overhang leftOverhang, rightOverhang ; // leftOverhangClassifier is for the overhang subexon at the left side of this intron.
} ;

struct _subexon
{
	char chrId ;
	int start, end ;
	int leftType, rightType ;
	double avgDepth ;
	//double ratio, classifier ;
	double leftRatio, rightRatio ;
	double leftClassifier, rightClassifier ;
	int lcCnt, rcCnt ;
	int leftStrand, rightStrand ;
	
	int nextCnt, prevCnt ;
	int *next, *prev ;

	// The information for the intron followed by this subexons, some of the
	// subexons may not followed by intron so we use pointer to save memory.
} ;

struct _seInterval
{
	int chrId ;
	int start, end ;
	int type ; // 0-subexon, 1-intronicInfo
	int idx ;
} ;

struct _subexonSplit
{
	char chrId ;
	int pos ;
	int type ; //1-start of a subexon. 2-end of a subexon 
	int splitType ; //0-soft boundary, 1-start of an exon, 2-end of an exon.
} ;
char buffer[4096] ;

bool CompSubexonSplit( struct _subexonSplit a, struct _subexonSplit b )
{
	if ( a.chrId < b.chrId )
		return true ;
	else if ( a.chrId > b.chrId )
		return false ;
	else if ( a.pos != b.pos )
		return a.pos < b.pos ;
	else
		return a.type < b.type ;
}

int InputSubexon( char *in, Alignments &alignments, struct _subexon &se, bool needPrevNext = false )
{
	int i, k ;
	char chrName[50] ;
	char ls[3], rs[3] ;	
	sscanf( in, "%s %d %d %d %d %s %s %lf %lf %lf %lf %lf", chrName, &se.start, &se.end, &se.leftType, &se.rightType, ls, rs,
			&se.avgDepth, &se.leftRatio, &se.rightRatio, 
			&se.leftClassifier, &se.rightClassifier ) ;	
	se.chrId = alignments.GetChromIdFromName( chrName ) ;
	se.nextCnt = se.prevCnt = 0 ;
	se.next = se.prev = NULL ;

	if ( ls[0] == '+' )
		se.leftStrand = 1 ;
	else if ( ls[0] == '-' )
		se.leftStrand = -1 ;
	else
		se.leftStrand = 0 ;
	
	if ( rs[0] == '+' )
		se.rightStrand = 1 ;
	else if ( rs[0] == '-' )
		se.rightStrand = -1 ;
	else
		se.rightStrand = 0 ;

	if ( needPrevNext )
	{
		char *p = in ;
		// Locate the offset for prevCnt
		for ( i = 0 ; i < 10 ; ++i )
		{
			p = strchr( p, ' ' ) ;
		}
		++p ;

		sscanf( p, "%d", &se.prevCnt ) ;
		p = strchr( p, ' ' ) ;
		++p ;
		se.prev = new int[ se.prevCnt ] ;
		for ( i = 0 ; i < se.prevCnt ; ++i )
		{
			sscanf( p, "%d", &se.prev[i] ) ;
			p = strchr( p, ' ' ) ;
			++p ;
		}

		sscanf( p, "%d", &se.nextCnt ) ;
		p = strchr( p, ' ' ) ;
		++p ;
		se.next = new int[ se.nextCnt ] ;
		for ( i = 0 ; i < se.nextCnt ; ++i )
		{
			sscanf( p, "%d", &se.next[i] ) ;
			p = strchr( p, ' ' ) ;
			++p ;
		}
	}
	return 1 ;
}

double GetUpdateMixtureGammeClassifier( double ratio, double cov, double piRatio, double kRatio[2], double thetaRatio[2],
	double piCov, double kCov[2], double thetaCov[2] )
{
	double p1 = 0, p2 ;
	if ( ratio > 0 )
		p1 = MixtureGammaAssignment( ratio, piRatio, kRatio, thetaRatio ) ;
	
	p2 = MixtureGammaAssignment( cov - 1, piCov, kCov, thetaCov ) ;
	double ret = 0 ;
	if ( p1 >= p2 ) // we should use ratio.
		ret = LogGammaDensity( ratio, kRatio[1], thetaRatio[1] ) 
			- LogGammaDensity( ratio, kRatio[0], thetaRatio[0] ) ;
	else
		ret = LogGammaDensity( cov - 1, kCov[1], thetaCov[1] ) 	
			- LogGammaDensity( cov - 1, kCov[0], thetaCov[0] ) ;

}

char StrandNumToSymbol( int strand )
{
	if ( strand > 0 )
		return '+' ;
	else if ( strand < 0 )
		return '-' ;
	else
		return '.' ;
}

int main( int argc, char *argv[] )
{
	int i, j, k ;
	FILE *fp ;
	std::vector<char *> files ;
	Alignments alignments ;

	if ( argc == 1 )
	{
		printf( "%s", usage ) ;
		return 0 ;
	}

	for ( i = 1 ; i < argc ; ++i )
	{
		if ( !strcmp( argv[i], "-s" ) )
		{
			files.push_back( argv[i + 1] ) ;
			++i ;
			continue ;
		}
	}
	int fileCnt = files.size() ;
	// Obtain the chromosome ids through bam file.
	fp = fopen( files[0], "r" ) ;		
	if ( fgets( buffer, sizeof( buffer ), fp ) != NULL )
	{
		int len = strlen( buffer ) ;
		buffer[len - 1] = '\0' ;
		alignments.Open( buffer + 1 ) ;
	}
	fclose( fp ) ;

	// Collect the split sites of subexons.
	std::vector<struct _subexonSplit> subexonSplits ;
	for ( k = 0 ; k < fileCnt ; ++k )
	{
		fp = fopen( files[k], "r" ) ;		
		struct _subexon se ;
		struct _subexonSplit sp ;
		char chrName[50] ;
		while ( fgets( buffer, sizeof( buffer), fp  ) != NULL )
		{
			if ( buffer[0] == '#' )
				continue ;

			InputSubexon( buffer, alignments, se ) ;
			// Ignore overhang subexons and ir subexons for now.
			if ( ( se.leftType == 0 && se.rightType == 1 ) 
				|| ( se.leftType == 2 && se.rightType == 0 ) 
				|| ( se.leftType == 2 && se.rightType == 1 ) )
				continue ;

			sp.chrId = se.chrId ;
			sp.pos = se.start ;
			sp.type = 1 ;
			sp.splitType = se.leftType ;				
			subexonSplits.push_back( sp ) ;

			sp.chrId = se.chrId ;
			sp.pos = se.end ;
			sp.type = 2 ;
			sp.splitType = se.rightType ;				
			subexonSplits.push_back( sp ) ;
		}	
		fclose( fp ) ;
	}
	// Pair up the split sites to get subexons
	std::sort( subexonSplits.begin(), subexonSplits.end(), CompSubexonSplit ) ;
	int splitCnt = subexonSplits.size() ;
	std::vector<struct _subexon> subexons ;
	int diffCnt = 0 ; // |start of subexon split| - |end of subexon split|

	for ( i = 0 ; i < splitCnt - 1 ; ++i )	
	{
		struct _subexon se ;
		/*if ( subexonSplits[i].pos == 553989 )
		{
			printf( "%d %d: %d %d\n", subexonSplits[i].pos, subexonSplits[i].type, subexonSplits[i + 1].pos, subexonSplits[i + 1].type ) ;
		}*/

		if ( subexonSplits[i].type == 1 )
			++diffCnt ;
		else
			--diffCnt ;

		if ( subexonSplits[i + 1].chrId != subexonSplits[i].chrId )
		{
			diffCnt = 0 ;		
			continue ;
		}

		if ( diffCnt == 0 ) // the interval between subexon
			continue ;

		se.chrId = subexonSplits[i].chrId ;
		se.start = subexonSplits[i].pos ;
		se.leftType = subexonSplits[i].splitType ;
		if ( subexonSplits[i].type == 2 )	
			++se.start ;

		se.end = subexonSplits[i + 1].pos ;
		se.rightType = subexonSplits[i + 1].splitType ;
		if ( subexonSplits[i + 1].type == 1 )
			--se.end ;
			
		if ( se.start > se.end ) //Note: this handles the case of repeated subexon split.
			continue ;
		se.leftClassifier = se.rightClassifier = 0 ;
		se.lcCnt = se.rcCnt = 0 ;
		subexons.push_back( se ) ;
	}
	// Merge the adjacent soft boundaries 
	std::vector<struct _subexon> rawSubexons = subexons ;
	int seCnt = subexons.size() ;
	subexons.clear() ;
	for ( i = 1, k = 0 ; i < seCnt ; ++i )
	{
		if ( rawSubexons[k].rightType == 0 && rawSubexons[i].leftType == 0 
			&& rawSubexons[k].end + 1 == rawSubexons[i].start )			
		{
			rawSubexons[k].end = rawSubexons[i].end ;
			rawSubexons[k].rightType = rawSubexons[k].rightType ;
		}
		else
		{
			subexons.push_back( rawSubexons[k] ) ;		
			k = i ;
		}
	}
	subexons.push_back( rawSubexons[k] ) ;		
	
	// Create the dummy intervals.
	seCnt = subexons.size() ;
	std::vector<struct _intronicInfo> intronicInfos ;
	std::vector<struct _seInterval> seIntervals ;
	for ( i = 0 ; i < seCnt - 1 ; ++i )
	{
		struct _seInterval ni ; // new interval
		ni.start = subexons[i].start ;
		ni.end = subexons[i].end ;
		ni.type = 0 ;
		ni.idx = i ;
		ni.chrId = subexons[i].chrId ;
		seIntervals.push_back( ni ) ;
		if ( subexons[i].chrId == subexons[i + 1].chrId && 
			subexons[i].end + 1 < subexons[i + 1].start &&
			subexons[i].rightType + subexons[i].leftType != 0 )
		{
			// Only consider the intervals like ]..[,]...(, )...[
			// The case like ]...] is actaully things like ][...] in subexon perspective,
			// so they won't pass the if-statement
			struct _intronicInfo nii ; // new intronic info
			ni.start = subexons[i].end + 1 ;
			ni.end = subexons[i + 1].start - 1 ;
			ni.type = 1 ;
			ni.idx = intronicInfos.size() ;
			
			nii.chrId = subexons[i].chrId ;
			nii.start = ni.start ;
			nii.end = ni.end ; 
			nii.irClassifier = 0 ;
			nii.irCnt = 0 ;
			nii.leftOverhang.cnt = 0 ;
			nii.leftOverhang.length = 0 ;
			nii.leftOverhang.classifier = 0 ;
			nii.rightOverhang.cnt = 0 ;
			nii.rightOverhang.length = 0 ;
			nii.rightOverhang.classifier = 0 ;
			intronicInfos.push_back( nii ) ;
		}
	}

	// Go through all the files to get some statistics number
	double avgIrPiRatio = 0 ;
	double avgIrPiCov = 0 ;
	double irPiRatio, irKRatio[2], irThetaRatio[2] ; // Some statistical results
	double irPiCov, irKCov[2], irThetaCov[2] ;
	
	double avgOverhangPiRatio = 0 ;
	double avgOverhangPiCov = 0 ;
	double overhangPiRatio, overhangKRatio[2], overhangThetaRatio[2] ; // Some statistical results
	double overhangPiCov, overhangKCov[2], overhangThetaCov[2] ;

	for ( k = 0 ; k < fileCnt ; ++k )
	{
		fp = fopen( files[k], "r" ) ;		
		
		while ( fgets( buffer, sizeof( buffer), fp  ) != NULL )
		{
			if ( buffer[0] == '#' )
			{
				char buffer2[100] ;
				sscanf( buffer, "%s", buffer2 ) ;	
				if ( !strcmp( buffer2, "#fitted_ir_parameter_ratio:" ) )
				{
					// TODO: ignore certain samples if the coverage seems wrong.
					sscanf( buffer, "%s %s %lf %s %lf %s %lf %s %lf %s %lf", 
								buffer2, buffer2, &irPiRatio, buffer2, &irKRatio[0], buffer2, &irThetaRatio[0],
								buffer2, &irKRatio[1], buffer2, &irThetaRatio[1] ) ;	
					avgIrPiRatio += irPiRatio ;
				}
				else if ( !strcmp( buffer2, "#fitted_ir_parameter_cov:" ) )
				{
				}
				else if ( !strcmp( buffer2, "#fitted_overhang_parameter_ratio:" ) )
				{
					sscanf( buffer, "%s %s %lf %s %lf %s %lf %s %lf %s %lf", 
								buffer2, buffer2, &overhangPiRatio, buffer2, &overhangKRatio[0], buffer2, &overhangThetaRatio[0],
								buffer2, &overhangKRatio[1], buffer2, &overhangThetaRatio[1] ) ;	
					avgOverhangPiRatio += overhangPiRatio ;
				}
			}
			else
				break ;
		}
	}
	avgIrPiRatio /= fileCnt ;
	avgOverhangPiRatio /= fileCnt ;

	// Go through all the files to put statistical results into each subexon.
	std::vector< struct _subexon > sampleSubexons ;
	int subexonCnt = subexons.size() ;
	for ( k = 0 ; k < fileCnt ; ++k )
	{
		fp = fopen( files[k], "r" ) ;		
		struct _subexon se ;
		struct _subexonSplit sp ;
		char chrName[50] ;
			
		sampleSubexons.clear() ;

		int tag = 0 ;
		while ( fgets( buffer, sizeof( buffer), fp  ) != NULL )
		{
			if ( buffer[0] == '#' )
			{
				char buffer2[200] ;
				sscanf( buffer, "%s", buffer2 ) ;	
				if ( !strcmp( buffer2, "#fitted_ir_parameter_ratio:" ) )
				{
					sscanf( buffer, "%s %s %lf %s %lf %s %lf %s %lf %s %lf", 
								buffer2, buffer2, &irPiRatio, buffer2, &irKRatio[0], buffer2, &irThetaRatio[0],
								buffer2, &irKRatio[1], buffer2, &irThetaRatio[1] ) ;	
				}
				else if ( !strcmp( buffer2, "#fitted_ir_parameter_cov:" ) )
				{
					sscanf( buffer, "%s %s %lf %s %lf %s %lf %s %lf %s %lf", 
								buffer2, buffer2, &irPiCov, buffer2, &irKCov[0], buffer2, &irThetaCov[0],
								buffer2, &irKCov[1], buffer2, &irThetaCov[1] ) ;	
				}
				else if ( !strcmp( buffer2, "#fitted_overhang_parameter_ratio:" ) )
				{
					sscanf( buffer, "%s %s %lf %s %lf %s %lf %s %lf %s %lf", 
								buffer2, buffer2, &overhangPiRatio, buffer2, &overhangKRatio[0], buffer2, &overhangThetaRatio[0],
								buffer2, &overhangKRatio[1], buffer2, &overhangThetaRatio[1] ) ;	
				}	
				else if ( !strcmp( buffer2, "#fitted_overhang_parameter_cov:" ) )
				{
					sscanf( buffer, "%s %s %lf %s %lf %s %lf %s %lf %s %lf", 
								buffer2, buffer2, &overhangPiCov, buffer2, &overhangKCov[0], buffer2, &overhangThetaCov[0],
								buffer2, &overhangKCov[1], buffer2, &overhangThetaCov[1] ) ;	
				}
				continue ;
			}

			InputSubexon( buffer, alignments, se ) ;
			sampleSubexons.push_back( se ) ;

			// TODO: add connection information.
		}
		
		int sampleSubexonCnt = sampleSubexons.size() ;
		int intervalCnt = seIntervals.size() ;
		for ( i = 0 ; i < sampleSubexonCnt ; ++i )	
		{
			struct _subexon &se = sampleSubexons[i] ;
			while ( tag < intervalCnt )	
			{
				if ( subexons[tag].chrId < se.chrId || subexons[tag].end < se.start )
				{
					++tag ;
					continue ;
				}
				else
					break ;
			}
			
			for ( j = tag ; j < intervalCnt ; ++j )
			{
				if ( seIntervals[j].start > se.end || seIntervals[j].chrId > se.chrId )
					break ;
				int idx ;	
				if ( seIntervals[j].type == 0 )
				{
					idx = seIntervals[j].idx ;
					if ( subexons[idx].leftType == 1 && se.leftType == 1 && subexons[idx].start == se.start )
					{
						subexons[idx].leftClassifier -= 2.0 * log( se.leftClassifier ) ;		
						++subexons[idx].lcCnt ;
					}
					if ( subexons[idx].rightType == 2 && se.rightType == 2 && subexons[idx].end == se.end )
					{
						subexons[idx].rightClassifier -= 2.0 * log( se.rightClassifier ) ;
						++subexons[idx].rcCnt ;
					}
				}
				else if ( seIntervals[j].type == 1 )
				{
					idx = seIntervals[j].idx ;
					// Overlap on the left part of intron
					if ( se.start < intronicInfos[idx].start && se.end < intronicInfos[idx].end )
					{
						int len = se.end - intronicInfos[idx].start + 1 ;
						intronicInfos[idx].leftOverhang.length += len ;
						++intronicInfos[idx].leftOverhang.cnt ;
						
						// Note that the sample subexon must have a soft boundary at right hand side, 
						// otherwise, this won't show up in intronic Info
						double update = GetUpdateMixtureGammeClassifier( se.leftRatio, se.avgDepth, 
							overhangPiRatio, overhangKRatio, overhangThetaRatio, 
							overhangPiCov, overhangKCov, overhangThetaCov ) ;
						intronicInfos[idx].leftOverhang.classifier += update ;				
					}
					// Overlap on the right part of intron
					else if ( se.start > intronicInfos[idx].start && se.end > intronicInfos[idx].end )
					{
						int len = intronicInfos[idx].end - se.start + 1 ;
						intronicInfos[idx].rightOverhang.length += len ;
						++intronicInfos[idx].rightOverhang.cnt ;
						
						// Note that the sample subexon must have a soft boundary at left hand side, 
						// otherwise, this won't show up in intronic Info
						double update = GetUpdateMixtureGammeClassifier( se.rightRatio, se.avgDepth, 
							overhangPiRatio, overhangKRatio, overhangThetaRatio, 
							overhangPiCov, overhangKCov, overhangThetaCov ) ;
						intronicInfos[idx].rightOverhang.classifier += update ;				
					}
					// Intron is fully contained in this sample subexon, then it is a ir candidate
					else if ( se.start <= intronicInfos[idx].start && se.end >= intronicInfos[idx].end )
					{
						if ( se.leftType == 2 && se.rightType == 1 )		
						{
							double update = GetUpdateMixtureGammeClassifier( se.leftRatio, se.avgDepth,
								irPiRatio, irKRatio, irThetaRatio,
								irPiCov, irKCov, irThetaCov ) ;
							intronicInfos[idx].irClassifier += update ;
							++intronicInfos[idx].irCnt ;
						}
						else if ( se.leftType == 1 && se.rightType == 2 )
						{
							intronicInfos[idx].irClassifier += LogGammaDensity( 4.0, irKRatio[1], irThetaRatio[1] )
							                                         - LogGammaDensity( 4.0, irKRatio[0], irThetaRatio[0] ) ;
							++intronicInfos[idx].irCnt ;
						}
						else
						{
							// the intron is contained in a overhang subexon from the sample or single-exon island
						}
					}
					// sample subexon is contained in the intron.
					else
					{
						// Do nothing.			
					}
				}
			}
		}
		fclose( fp ) ;
	}

	// Convert the temporary statistics number into formal statistics result.
	for ( i = 0 ; i < subexonCnt ; ++i ) 
	{
		struct _subexon &se = subexons[i] ;
		se.leftClassifier = 1 - chicdf( se.leftClassifier, 2 * se.lcCnt ) ; 	
		se.rightClassifier = 1 - chicdf( se.rightClassifier, 2 * se.rcCnt ) ;
	}

	int iiCnt = intronicInfos.size() ; //intronicInfo count
	for ( i = 0 ; i < iiCnt ; ++i )
	{
		struct _intronicInfo &ii = intronicInfos[i] ;
		ii.irClassifier = (double)1.0 / ( 1.0 + exp( ii.irClassifier + log( 1 - avgIrPiRatio ) - log( avgIrPiRatio ) ) ) ;
		ii.leftOverhang.classifier = (double)1.0 / ( 1.0 + exp( ii.leftOverhang.classifier + 
					log( 1 - avgOverhangPiRatio ) - log( avgOverhangPiRatio ) ) ) ;
		ii.rightOverhang.classifier = (double)1.0 / ( 1.0 + exp( ii.rightOverhang.classifier + 
					log( 1 - avgOverhangPiRatio ) - log( avgOverhangPiRatio ) ) ) ;
	}

	// Change the classifier for the hard boundaries if its adjacent intron has intron retention classifier
	int intervalCnt = seIntervals.size() ;
	for ( i = 0 ; i < intervalCnt ; ++i )
	{
		if ( seIntervals[i].type == 1 && intronicInfos[ seIntervals[i].idx ].irCnt > 0 )
		{
			int idx = seIntervals[i].idx ;
			if ( intronicInfos[idx].leftOverhang.cnt > 0 )
			{
				int k = seIntervals[i - 1].idx ;
				// Should aim for more conservative?
				if ( subexons[k].rightClassifier > intronicInfos[idx].leftOverhang.classifier )
					subexons[k].rightClassifier = intronicInfos[idx].leftOverhang.classifier ;
			}

			if ( intronicInfos[idx].rightOverhang.cnt > 0 )
			{
				int k = seIntervals[i + 1].idx ;
				if ( subexons[k].leftClassifier > intronicInfos[idx].rightOverhang.classifier )
					subexons[k].leftClassifier = intronicInfos[idx].rightOverhang.classifier ;
			}
		}
	}

	// Output the result.
	for ( i = 0 ; i < intervalCnt ; ++i )
	{
		if ( seIntervals[i].type == 0 )
		{
			struct _subexon &se = subexons[i] ;
			
			char ls, rs ;
			if ( se.leftStrand == 1 )
				ls = '+' ;
			else if ( se.leftStrand == -1 )
				ls = '-' ;
			else
				ls = '.' ;

			if ( se.rightStrand == 1 )
				rs = '+' ;
			else if ( se.rightStrand == -1 )
				rs = '-' ;
			else
				rs = '.' ;

			printf( "%s %d %d %d %d %c %c -1 -1 -1 %lf %lf ", alignments.GetChromName( se.chrId ), se.start, se.end,
					se.leftType, se.rightType, ls, rs, se.leftClassifier, se.rightClassifier ) ;
			if ( i > 0 && seIntervals[i - 1].chrId == seIntervals[i].chrId 
				&& seIntervals[i - 1].end + 1 == seIntervals[i].start 
				&& !( seIntervals[i - 1].type == 0 && 
					subexons[ seIntervals[i - 1].idx ].rightType != se.leftType ) 
				&& !( seIntervals[i - 1].type == 1 && intronicInfos[ seIntervals[i - 1].idx ].irCnt == 0
					&& intronicInfos[ seIntervals[i - 1].idx ].rightOverhang.cnt == 0 ) )
			{
				printf( "%d ", se.prevCnt + 1 ) ;
				for ( j = 0 ; j < se.prevCnt ; ++j )
					printf( "%d ", se.prev[j] ) ;
				printf( "%d ", se.start - 1 ) ;
			}
			else
			{
				printf( "%d ", se.prevCnt ) ;
				for ( j = 0 ; j < se.prevCnt ; ++j )
					printf( "%d ", se.prev[j] ) ;
			}

			if ( i < intervalCnt - 1 && seIntervals[i].chrId == seIntervals[i + 1].chrId 
				&& seIntervals[i].end == seIntervals[i + 1].start - 1
				&& !( seIntervals[i + 1].type == 0 &&
					subexons[ seIntervals[i + 1].idx ].leftType != se.rightType ) 
				&& !( seIntervals[i + 1].type == 1 && intronicInfos[ seIntervals[i + 1].idx ].irCnt == 0
					&& intronicInfos[ seIntervals[i + 1].idx ].leftOverhang.cnt == 0 ) )
			{
				printf( "%d %d ", se.nextCnt + 1, se.end + 1 ) ;
			}
			else
				printf( "%d ", se.nextCnt ) ;
			for ( j = 0 ; j < se.nextCnt ; ++j )
				printf( "%d ", se.next[j] ) ;
			printf( "\n" ) ;
		}
		else if ( seIntervals[i].type == 1 )
		{
			struct _intronicInfo &ii = intronicInfos[ seIntervals[i].idx ] ;
			if ( ii.irCnt > 0 )
			{
				printf( "%s %d %d 2 1 %c %c -1 -1 -1 %lf %lf 1 %d 1 %d\n",
					alignments.GetChromName( ii.chrId ), ii.start, ii.end, 
					StrandNumToSymbol( subexons[ seIntervals[i - 1].idx ].rightStrand ), 
					StrandNumToSymbol( subexons[ seIntervals[i + 1].idx ].leftStrand ), 
					ii.irClassifier, ii.irClassifier,
					seIntervals[i - 1].end, seIntervals[i + 1].start ) ;
			}
			else
			{
				// left overhang.
				if ( ii.leftOverhang.cnt > 0 )
				{
					printf( "%s %d %d 2 0 %c . -1 -1 -1 %lf %lf 1 %d 0\n",
						alignments.GetChromName( ii.chrId ), ii.start, 
						ii.start + ( ii.leftOverhang.length /  ii.leftOverhang.cnt ) - 1,
						StrandNumToSymbol( subexons[ seIntervals[i - 1].idx ].rightStrand ),
						ii.leftOverhang.classifier, ii.leftOverhang.classifier,
						ii.start - 1 ) ; 
				}

				// right overhang.
				if ( ii.rightOverhang.cnt > 0 )
				{
					printf( "%s %d %d 0 1 . %c -1 -1 -1 %lf %lf 0 0 %d\n",
						alignments.GetChromName( ii.chrId ), 
						ii.end - ( ii.rightOverhang.length / ii.rightOverhang.cnt ) + 1, ii.end,
						StrandNumToSymbol( subexons[ seIntervals[i + 1 ].idx ].leftStrand ),
						ii.rightOverhang.classifier, ii.rightOverhang.classifier,
						ii.end + 1 ) ;
				}

			}
		}
	}

	return 0 ;
}


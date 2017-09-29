#include <stdio.h>
#include <string.h>

#include <algorithm>
#include <vector>
#include <map>
#include <string>

#include "alignments.hpp"
#include "blocks.hpp"
#include "stats.hpp"
#include "SubexonGraph.hpp"

char usage[] = "combineSubexons [options]\n"
	       "Required options:\n"
	       "\t-s STRING: the path to the predicted subexon information. Can use multiple -s to specify multiple subexon prediction files\n" 
	       "\t\tor\n"
	       "\t--ls STRING: the path to file of the list of the predicted subexon information.\n" ;

struct _overhang
{
	int cnt ; //  the number of samples support this subexon.
	int validCnt ; // The number of samples that are used for compute probability.
	int length ;
	double classifier ;
} ;

struct _intronicInfo
{
	int chrId ;
	int start, end ;
	double irClassifier ;
	int irCnt ;
	int validIrCnt ;
	struct _overhang leftOverhang, rightOverhang ; // leftOverhangClassifier is for the overhang subexon at the left side of this intron.
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
	int chrId ;
	int pos ;
	int type ; //1-start of a subexon. 2-end of a subexon 
	int splitType ; //0-soft boundary, 1-start of an exon, 2-end of an exon.
	int strand ;
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
	else if ( a.type != b.type )
	{
		// the split site with no strand information should come first.
		if ( a.type == 0 ) 
			return true ;
		else if ( b.type == 0 )
			return false ;
		return a.type < b.type ;
	}
	else if ( a.splitType != b.splitType )
	{
		return a.splitType < b.splitType ;
	}
	
	return false ;
}

bool CompIrFromSamples( struct _seInterval a, struct _seInterval b )
{
	if ( a.chrId < b.chrId )
		return true ;
	else if ( a.chrId > b.chrId )
		return false ;
	else if ( a.start < b.start )
		return true ;
	else if ( a.start > b.start )
		return false ;
	else if ( a.end < b.end )
		return true ;
	else
		return false ;
}

// Keep this the same as in SubexonInfo.cpp.
double TransformCov( double c )
{
	return sqrt( c ) - 1 ;
}

double GetUpdateMixtureGammeClassifier( double ratio, double cov, double piRatio, double kRatio[2], double thetaRatio[2],
	double piCov, double kCov[2], double thetaCov[2] )
{
	double p1 = 0, p2 ;

	cov = TransformCov( cov ) ;
	if ( cov < ( kCov[0] - 1 ) * thetaCov[0] )
		cov = ( kCov[0] - 1 ) * thetaCov[0] ;

	if ( ratio > 0 )
		p1 = MixtureGammaAssignment( ratio, piRatio, kRatio, thetaRatio ) ;
	// Make sure cov > 1?	
	p2 = MixtureGammaAssignment( cov, piCov, kCov, thetaCov ) ;
	double ret = 0 ;
	if ( p1 >= p2 ) // we should use ratio.
		ret = LogGammaDensity( ratio, kRatio[1], thetaRatio[1] ) 
			- LogGammaDensity( ratio, kRatio[0], thetaRatio[0] ) ;
	else
		ret = LogGammaDensity( cov, kCov[1], thetaCov[1] ) 	
			- LogGammaDensity( cov, kCov[0], thetaCov[0] ) ;
	return ret ;
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

int StrandSymbolToNum( char c )
{
	if ( c == '+' )
		return 1 ;
	else if ( c == '-' )
		return -1 ;
	else
		return 0 ;
}

int *MergePositions( int *old, int ocnt, int *add, int acnt, int &newCnt )
{
	int i, j, k ;
	int *ret ;
	if ( acnt == 0 )	
	{
		newCnt = ocnt ;
		return old ;
	}
	if ( ocnt == 0 )
	{
		newCnt = acnt ;
		ret = new int[acnt] ;
		for ( i = 0 ; i < acnt ; ++i )
		{
			ret[i] = add[i] ;		
		}
		return ret ;
	}
	newCnt = 0 ;
	for ( i = 0, j = 0 ; i < ocnt && j < acnt ; )
	{
		if ( old[i] < add[j] )
		{
			++i ;
			++newCnt ;
		}
		else if ( old[i] == add[j] )
		{
			++i ; ++j ;
			++newCnt ;
		}
		else 
		{
			++j ;
			++newCnt ;
		}
	}
	newCnt = newCnt + ( ocnt - i ) + ( acnt - j ) ;
	// no new elements.
	if ( newCnt == ocnt )
		return old ;
	k = 0 ;
	//delete []old ;
	ret = new int[ newCnt ] ;
	for ( i = 0, j = 0 ; i < ocnt && j < acnt ; )
	{
		if ( old[i] < add[j] )
		{
			ret[k] = old[i] ;
			++i ;
			++k ;
		}
		else if ( old[i] == add[j] )
		{
			ret[k] = old[i] ;
			++i ; ++j ;
			++k ;
		}
		else 
		{
			ret[k] = add[j] ;
			++j ;
			++k ;
		}
	}
	for ( ; i < ocnt ; ++i, ++k )
		ret[k] = old[i] ;
	for ( ; j < acnt ; ++j, ++k )
		ret[k] = add[j] ;
	return ret ;
}



int main( int argc, char *argv[] )
{
	int i, j, k ;
	FILE *fp ;
	std::vector<char *> files ;

	Blocks regions ;
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
		else if ( !strcmp( argv[i], "--ls" ) )
		{
			FILE *fpLs = fopen( argv[i + 1], "r" ) ;
			char buffer[1024] ;
			while ( fgets( buffer, sizeof( buffer ), fpLs ) != NULL )
			{
				int len = strlen( buffer ) ;
				if ( buffer[len - 1] == '\n' )
				{
					buffer[len - 1] = '\0' ;
					--len ;

				}
				char *fileName = strdup( buffer ) ;
				files.push_back( fileName ) ;
			}
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
	std::vector<struct _seInterval> irFromSamples ;

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

			SubexonGraph::InputSubexon( buffer, alignments, se ) ;

			// Record all the intron rentention from the samples
			if ( se.leftType == 2 && se.rightType == 1 )
			{
				struct _seInterval si ;
				si.chrId = se.chrId ;
				si.start = se.start ;
				si.end = se.end ;

				irFromSamples.push_back( si ) ;
			}
			
			// Ignore overhang subexons and ir subexons for now.
			if ( ( se.leftType == 0 && se.rightType == 1 ) 
				|| ( se.leftType == 2 && se.rightType == 0 ) 
				|| ( se.leftType == 2 && se.rightType == 1 ) )
				continue ;
			
			sp.chrId = se.chrId ;
			sp.pos = se.start ;
			sp.type = 1 ;
			sp.splitType = se.leftType ;			
			sp.strand = se.leftStrand ;
			subexonSplits.push_back( sp ) ;

			sp.chrId = se.chrId ;
			sp.pos = se.end ;
			sp.type = 2 ;
			sp.splitType = se.rightType ;				
			sp.strand = se.rightStrand ;
			subexonSplits.push_back( sp ) ;
		}	
		fclose( fp ) ;
	}
	// Pair up the split sites to get subexons
	std::sort( subexonSplits.begin(), subexonSplits.end(), CompSubexonSplit ) ;

	// Force the soft boundary that collides with hard boundaries to be hard boundary.
	int splitCnt = subexonSplits.size() ;
	for ( i = 0 ; i < splitCnt ; ++i )
	{
		if ( subexonSplits[i].splitType != 0 )
			continue ;
		int newSplitType = 0 ;
		int newStrand = subexonSplits[i].strand ;
		for ( j = i + 1 ; j < splitCnt ; ++j )
		{
			if ( subexonSplits[i].type != subexonSplits[j].type || subexonSplits[i].pos != subexonSplits[j].pos ||
					subexonSplits[i].chrId != subexonSplits[j].chrId )
				break ;
			if ( subexonSplits[j].splitType != 0 )
			{
				newSplitType = subexonSplits[j].splitType ;
				newStrand = subexonSplits[j].strand ;
				break ;
			}
		}
		subexonSplits[i].splitType = newSplitType ;
		subexonSplits[i].strand = newStrand ;
	}
	
	// Build subexons from the collected split sites.
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
		se.leftStrand = subexonSplits[i].strand ;
		if ( subexonSplits[i].type == 2 )	
			++se.start ;

		se.end = subexonSplits[i + 1].pos ;
		se.rightType = subexonSplits[i + 1].splitType ;
		se.rightStrand = subexonSplits[i + 1].strand ;
		if ( subexonSplits[i + 1].type == 1 )
			--se.end ;
			
		if ( se.start > se.end ) //Note: this handles the case of repeated subexon split.
			continue ;
		se.leftClassifier = se.rightClassifier = 0 ;
		se.lcCnt = se.rcCnt = 0 ;
	
		se.next = se.prev = NULL ;
		se.nextCnt = se.prevCnt = 0 ;
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
			rawSubexons[k].rightType = rawSubexons[i].rightType ;
			rawSubexons[k].rightStrand = rawSubexons[i].rightStrand ;
		}
		else
		{
			subexons.push_back( rawSubexons[k] ) ;		
			k = i ;
		}
	}
	subexons.push_back( rawSubexons[k] ) ;		

	// Remove overhang, ir subexons intron created after putting multiple sample to gether.
	// eg: s0: [......)
	//     s1: [...]--------[....]
	//     s2: [...]..)-----[....]
	// Though the overhang from s2 is filtered in readin, there will a new overhang created combining s0,s1.
	// 	But be careful about how to compute the classifier for the overhang part contributed from s0.
	// Furthermore, note that the case of single-exon island showed up in intron retention region after combining is not possible when get here.
	//    eg: s0:[...]-----[...]
	//        s1:      (.)
	//        s2:[.............]
	//  After merge adjacent soft boundaries, the single-exon island will disappear.
	rawSubexons = subexons ;
	seCnt = subexons.size() ;
	subexons.clear() ;
	for ( i = 0 ; i < seCnt ; ++i )
	{
		if ( ( rawSubexons[i].leftType == 2 && rawSubexons[i].rightType == 1 )		// ir
			|| ( rawSubexons[i].leftType == 2 && rawSubexons[i].rightType == 0 )    // overhang	
			|| ( rawSubexons[i].leftType == 0 && rawSubexons[i].rightType == 1 ) )  
			continue ;
		subexons.push_back( rawSubexons[i] ) ;
	}
	
	// Remove the single-exon island if it shows up in the intron that is intron retentioned in some sample.
	rawSubexons = subexons ;
	seCnt = subexons.size() ;
	subexons.clear() ;
	k = 0 ;
	std::sort( irFromSamples.begin(), irFromSamples.end(), CompIrFromSamples ) ;
	int irFromSampleCnt = irFromSamples.size() ;

	for ( i = 0 ; i < seCnt ; ++i )
	{
		if ( rawSubexons[i].leftType != 0 || rawSubexons[i].rightType != 0 )
		{
			subexons.push_back( rawSubexons[i] ) ;
			continue ;
		}
		
		while ( k < irFromSampleCnt )
		{
			// Locate the ir that ends after the island.
			if ( irFromSamples[k].chrId < rawSubexons[i].chrId 
				|| ( irFromSamples[k].chrId == rawSubexons[i].chrId && irFromSamples[k].end < rawSubexons[i].end ) )
			{
				++k ;
				continue ;
			}
			break ;
		}
		bool contained = false ;
		for ( j = k ; j < irFromSampleCnt ; ++j )
		{
			if ( irFromSamples[j].chrId > rawSubexons[i].chrId || irFromSamples[j].start > rawSubexons[i].start )
				break ;
			if ( irFromSamples[j].start <= rawSubexons[i].start || irFromSamples[j].end >= rawSubexons[i].end )
			{
				contained = true ;
				break ;
			}
		}

		if ( !contained )
			subexons.push_back( rawSubexons[i] ) ;
	}
	std::vector<struct _seInterval>().swap( irFromSamples ) ;
	
	// Create the dummy intervals.
	seCnt = subexons.size() ;
	std::vector<struct _intronicInfo> intronicInfos ;
	std::vector<struct _seInterval> seIntervals ;
	for ( i = 0 ; i < seCnt ; ++i )
	{
		struct _seInterval ni ; // new interval
		ni.start = subexons[i].start ;
		ni.end = subexons[i].end ;
		ni.type = 0 ;
		ni.idx = i ;
		ni.chrId = subexons[i].chrId ;
		seIntervals.push_back( ni ) ;
		
		/*int nexti ;
		for ( nexti = i + 1 ; nexti < seCnt ; ++nexti )
			if ( subexons[ nexti ].leftType == 0 && subexons[nexti].rightType == 0 )*/

		if ( i < seCnt - 1 && subexons[i].chrId == subexons[i + 1].chrId && 
			subexons[i].end + 1 < subexons[i + 1].start &&
			subexons[i].rightType + subexons[i + 1].leftType != 0 )
		{
			// Only consider the intervals like ]..[,]...(, )...[
			// The case like ]...] is actaully things like ][...] in subexon perspective,
			// so they won't pass the if-statement
			struct _intronicInfo nii ; // new intronic info
			ni.start = subexons[i].end + 1 ;
			ni.end = subexons[i + 1].start - 1 ;
			ni.type = 1 ;
			ni.idx = intronicInfos.size() ;
			seIntervals.push_back( ni ) ;
			
			nii.chrId = subexons[i].chrId ;
			nii.start = ni.start ;
			nii.end = ni.end ; 
			nii.irClassifier = 0 ;
			nii.irCnt = 0 ;
			nii.validIrCnt = 0 ;
			nii.leftOverhang.cnt = 0 ;
			nii.leftOverhang.validCnt = 0 ;
			nii.leftOverhang.length = 0 ;
			nii.leftOverhang.classifier = 0 ;
			nii.rightOverhang.cnt = 0 ;
			nii.rightOverhang.validCnt = 0 ;
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

			SubexonGraph::InputSubexon( buffer, alignments, se, true ) ;
			sampleSubexons.push_back( se ) ;
		}
		
		int sampleSubexonCnt = sampleSubexons.size() ;
		int intervalCnt = seIntervals.size() ;
		for ( i = 0 ; i < sampleSubexonCnt ; ++i )	
		{
			struct _subexon &se = sampleSubexons[i] ;
			while ( tag < intervalCnt )	
			{
				if ( seIntervals[tag].chrId < se.chrId || 
					( seIntervals[tag].chrId == se.chrId && seIntervals[tag].end < se.start ) )
				{
					++tag ;
					continue ;
				}
				else
					break ;
			}
			
			for ( j = tag ; j < intervalCnt ; ++j )
			{
				if ( seIntervals[j].start > se.end || seIntervals[j].chrId > se.chrId ) // terminate if no overlap.
					break ;
				int idx ;	
				
				if ( seIntervals[j].type == 0 )
				{
					idx = seIntervals[j].idx ;
					if ( subexons[idx].leftType == 1 && se.leftType == 1 && subexons[idx].start == se.start )
					{
						double tmp = se.leftClassifier ;
						if ( se.leftClassifier == 0 )
							tmp = 1e-7 ;
						subexons[idx].leftClassifier -= 2.0 * log( tmp ) ;		
						++subexons[idx].lcCnt ;

						subexons[idx].prev = MergePositions( subexons[idx].prev, subexons[idx].prevCnt, 
										se.prev, se.prevCnt, subexons[idx].prevCnt ) ;
					}
					if ( subexons[idx].rightType == 2 && se.rightType == 2 && subexons[idx].end == se.end )
					{
						double tmp = se.rightClassifier ;
						if ( se.rightClassifier == 0 )
							tmp = 1e-7 ;
						subexons[idx].rightClassifier -= 2.0 * log( tmp ) ;
						++subexons[idx].rcCnt ;

						subexons[idx].next = MergePositions( subexons[idx].next, subexons[idx].nextCnt, 
										se.next, se.nextCnt, subexons[idx].nextCnt ) ;
					}

					if ( subexons[idx].leftType == 0 && subexons[idx].rightType == 0
						&& se.leftType == 0 && se.rightType == 0 ) // the single-exon island.
					{
						double tmp = se.leftClassifier ;
						if ( se.leftClassifier == 0 )
							tmp = 1e-7 ;
						subexons[idx].leftClassifier -= 2.0 * log( tmp ) ;
						subexons[idx].rightClassifier = subexons[idx].leftClassifier ;
						++subexons[idx].lcCnt ;
						++subexons[idx].rcCnt ;
					}
				}
				else if ( seIntervals[j].type == 1 )
				{
					idx = seIntervals[j].idx ;
					// Overlap on the left part of intron
					if ( se.start <= intronicInfos[idx].start && se.end < intronicInfos[idx].end )
					{
						int len = se.end - intronicInfos[idx].start + 1 ;
						intronicInfos[idx].leftOverhang.length += len ;
						++intronicInfos[idx].leftOverhang.cnt ;
						
						// Note that the sample subexon must have a soft boundary at right hand side, 
						// otherwise, this part is not an intron and won't show up in intronic Info.
						if ( se.leftType == 2 )
						{
							if ( se.leftRatio > 0 && se.avgDepth > 1 )
							{
								++intronicInfos[idx].leftOverhang.validCnt ;

								double update = GetUpdateMixtureGammeClassifier( se.leftRatio, se.avgDepth, 
										overhangPiRatio, overhangKRatio, overhangThetaRatio, 
										overhangPiCov, overhangKCov, overhangThetaCov ) ;
								intronicInfos[idx].leftOverhang.classifier += update ;				
							}
						}
						else if ( se.leftType == 1 )
						{
							++intronicInfos[idx].leftOverhang.validCnt ;
							double update = GetUpdateMixtureGammeClassifier( 1.0, se.avgDepth, 
									overhangPiRatio, overhangKRatio, overhangThetaRatio, 
									overhangPiCov, overhangKCov, overhangThetaCov ) ;
							intronicInfos[idx].leftOverhang.classifier += update ;				

						}
						// ignore the contribution of single-exon island here?
					}
					// Overlap on the right part of intron
					else if ( se.start > intronicInfos[idx].start && se.end >= intronicInfos[idx].end )
					{
						int len = intronicInfos[idx].end - se.start + 1 ;
						intronicInfos[idx].rightOverhang.length += len ;
						++intronicInfos[idx].rightOverhang.cnt ;
						
						// Note that the sample subexon must have a soft boundary at left hand side, 
						// otherwise, this won't show up in intronic Info
						if ( se.rightType == 1 )
						{
							if ( se.rightRatio > 0 && se.avgDepth > 1 )
							{
								++intronicInfos[idx].rightOverhang.validCnt ;

								double update = GetUpdateMixtureGammeClassifier( se.rightRatio, se.avgDepth, 
										overhangPiRatio, overhangKRatio, overhangThetaRatio, 
										overhangPiCov, overhangKCov, overhangThetaCov ) ;
								intronicInfos[idx].rightOverhang.classifier += update ;				
							}
						}
						else if ( se.rightType == 2 )
						{
							++intronicInfos[idx].rightOverhang.validCnt ;

							double update = GetUpdateMixtureGammeClassifier( 1, se.avgDepth, 
									overhangPiRatio, overhangKRatio, overhangThetaRatio, 
									overhangPiCov, overhangKCov, overhangThetaCov ) ;
							intronicInfos[idx].rightOverhang.classifier += update ;				

						}
					}
					// Intron is fully contained in this sample subexon, then it is a ir candidate
					else if ( se.start <= intronicInfos[idx].start && se.end >= intronicInfos[idx].end )
					{
						if ( se.leftType == 2 && se.rightType == 1 )		
						{
							double ratio = regions.PickLeftAndRightRatio( se.leftRatio, se.rightRatio ) ;
							++intronicInfos[idx].irCnt ;
							if ( ratio > 0 && se.avgDepth > 1 )
							{
								double update = GetUpdateMixtureGammeClassifier( ratio, se.avgDepth,
										irPiRatio, irKRatio, irThetaRatio,
										irPiCov, irKCov, irThetaCov ) ;
								//if ( se.start == 17171 )
								//	printf( "hi %lf %d %d: %d %d\n", update, se.start, se.end, intronicInfos[idx].start, intronicInfos[idx].end ) ;
								intronicInfos[idx].irClassifier += update ;
								++intronicInfos[idx].validIrCnt ;
							}
						}
						else if ( se.leftType == 1 && se.rightType == 2 )
						{
							//intronicInfos[idx].irClassifier += LogGammaDensity( 4.0, irKRatio[1], irThetaRatio[1] )
							//                                         - LogGammaDensity( 4.0, irKRatio[0], irThetaRatio[0] ) ;
							/*if ( se.start == 6643539 )
							{
								printf( "%lf: %lf %lf\n", se.avgDepth, MixtureGammaAssignment( ( irKCov[0] - 1 ) * irThetaCov[0], irPiRatio, irKCov, irThetaCov ),
									MixtureGammaAssignment( TransformCov( 4.0 ), irPiRatio, irKCov, irThetaCov ) ) ;
							}*/
							if ( se.avgDepth > 1 )
							{
								// let the depth be the threshold to determine.
								double update = GetUpdateMixtureGammeClassifier( 4.0, se.avgDepth,
										irPiRatio, irKRatio, irThetaRatio,
										irPiCov, irKCov, irThetaCov ) ;
								intronicInfos[idx].irClassifier += update ;
								++intronicInfos[idx].irCnt ;
								++intronicInfos[idx].validIrCnt ;
							}
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
		
		for ( i = 0 ; i < sampleSubexonCnt ; ++i )
		{
			if ( sampleSubexons[i].nextCnt > 0 )
				delete[] sampleSubexons[i].next ;
			if ( sampleSubexons[i].prevCnt > 0 )
				delete[] sampleSubexons[i].prev ;
		}
	}
	
	// Convert the temporary statistics number into formal statistics result.
	for ( i = 0 ; i < subexonCnt ; ++i ) 
	{
		struct _subexon &se = subexons[i] ;
		if ( se.leftType == 0 && se.rightType == 0 ) // single-exon txpt.
		{
			se.leftClassifier = se.rightClassifier = 1 - chicdf( se.rightClassifier, 2 * se.rcCnt ) ;
		}
		else
		{
			if ( se.leftType == 1 )
				se.leftClassifier = 1 - chicdf( se.leftClassifier, 2 * se.lcCnt ) ; 	
			else
				se.leftClassifier = -1 ;

			if ( se.rightType == 2 )
				se.rightClassifier = 1 - chicdf( se.rightClassifier, 2 * se.rcCnt ) ;
			else
				se.rightClassifier = -1 ;
		}
	}

	int iiCnt = intronicInfos.size() ; //intronicInfo count
	for ( i = 0 ; i < iiCnt ; ++i )
	{
		struct _intronicInfo &ii = intronicInfos[i] ;
		if ( ii.validIrCnt > 0 )
			ii.irClassifier = (double)1.0 / ( 1.0 + exp( ii.irClassifier + log( 1 - avgIrPiRatio ) - log( avgIrPiRatio ) ) ) ;
		else
			ii.irClassifier = -1 ;
		
		if ( ii.leftOverhang.validCnt > 0 )
			ii.leftOverhang.classifier = (double)1.0 / ( 1.0 + exp( ii.leftOverhang.classifier + 
						log( 1 - avgOverhangPiRatio ) - log( avgOverhangPiRatio ) ) ) ;
		else
			ii.leftOverhang.classifier = -1 ;

		if ( ii.rightOverhang.validCnt > 0 )
			ii.rightOverhang.classifier = (double)1.0 / ( 1.0 + exp( ii.rightOverhang.classifier + 
						log( 1 - avgOverhangPiRatio ) - log( avgOverhangPiRatio ) ) ) ;
		else
			ii.rightOverhang.classifier = -1 ;
	}

	// Change the classifier for the hard boundaries if its adjacent intron has intron retention classifier
	//    which collide with overhang subexon.
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
			struct _subexon &se = subexons[ seIntervals[i].idx ] ;
			
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
					&& intronicInfos[ seIntervals[i - 1].idx ].rightOverhang.cnt == 0 ) 
				&& ( se.prevCnt == 0 || se.start - 1 != se.prev[ se.prevCnt - 1 ] ) ) // The connection showed up in the subexon file.
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
					&& intronicInfos[ seIntervals[i + 1].idx ].leftOverhang.cnt == 0 ) 
				&& ( se.nextCnt == 0 || se.end + 1 != se.next[0] ) )
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
				printf( "%s %d %d 2 1 . . -1 -1 -1 %lf %lf 1 %d 1 %d\n",
					alignments.GetChromName( ii.chrId ), ii.start, ii.end, 
					ii.irClassifier, ii.irClassifier,
					seIntervals[i - 1].end, seIntervals[i + 1].start ) ;
			}
			else
			{
				// left overhang.
				if ( ii.leftOverhang.cnt > 0 )
				{
					printf( "%s %d %d 2 0 . . -1 -1 -1 %lf %lf 1 %d 0\n",
						alignments.GetChromName( ii.chrId ), ii.start, 
						ii.start + ( ii.leftOverhang.length /  ii.leftOverhang.cnt ) - 1,
						ii.leftOverhang.classifier, ii.leftOverhang.classifier,
						ii.start - 1 ) ; 
				}

				// right overhang.
				if ( ii.rightOverhang.cnt > 0 )
				{
					printf( "%s %d %d 0 1 . . -1 -1 -1 %lf %lf 0 1 %d\n",
						alignments.GetChromName( ii.chrId ), 
						ii.end - ( ii.rightOverhang.length / ii.rightOverhang.cnt ) + 1, ii.end,
						ii.rightOverhang.classifier, ii.rightOverhang.classifier,
						ii.end + 1 ) ;
				}

			}
		}
	}

	return 0 ;
}


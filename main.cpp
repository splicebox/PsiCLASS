// report the total depth for the segment
// Format: ./a.out alignments.bam splice_site
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include <algorithm>
#include <vector>
#include <math.h>

#include "alignments.hpp"
#include "blocks.hpp"
#include "stats.hpp"

#define ABS(x) ((x)<0?-(x):(x))

char usage[] = "./subexon-info alignment.bam intron.splice [options]\n"
		"options:\n"
		"\t--minDepth INT: the minimum coverage depth considered as part of a subexon (default: 2)\n"
		"\t--noStats: do not compute the statistical scores (default: not used)\n" ;
char buffer[4096] ;

int gMinDepth ;

bool CompSplitSite( struct _splitSite a, struct _splitSite b )
{
	if ( a.chrId < b.chrId )
		return true ;
	else if ( a.chrId > b.chrId )
		return false ;
	else if ( a.pos != b.pos )
		return a.pos < b.pos ;
	else if ( a.type != b.type )
		return a.type < b.type ; // We want the start of exons comes first, since we are scanning the genome from left to right,
					 // so we can terminate the extension early and create a single-base subexon later.
	else
		return a.oppositePos < b.oppositePos ;
}

bool CompBlocksByAvgDepth( struct _block a, struct _block b )
{
	double avgA = a.depthSum / (double)( a.end - a.start + 1 ) ;
	double avgB = b.depthSum / (double)( b.end - b.start + 1 ) ;
	return avgA < avgB ;
}

bool CompBlocksByRatio( struct _block a, struct _block b )
{
	return a.ratio < b.ratio ;	
}

int CompDouble( const void *p1, const void *p2 )
{
	return *(double *)p1 - *(double *)p2 ;
}

void CleanAndSortSplitSites( std::vector< struct _splitSite> &sites )
{
	std::sort( sites.begin(), sites.end(), CompSplitSite ) ;	
	int i ;
	int size = sites.size() ;
	int k = 1 ;
	for ( i = 1 ; i < size ; ++i )
	{
		if ( sites[i].chrId != sites[k - 1].chrId || sites[i].pos != sites[k - 1].pos || sites[i].type != sites[k - 1].type )
		{
			sites[k] = sites[i] ;
			++k ;
		}
		/*else
		{
			if ( sites[i].type != sites[k-1].type )
			{
				printf( "%d\n", sites[i].pos ) ;
			}
		}*/
	}
	for ( i = size - 1 ; i >= k ; --i )
		sites.pop_back() ;

	// For the sites that corresponds to the start of an exon, we remove the adjust to it and 
	/*for ( i = 1 ; i < size ; ++i )
	{
		if ( sites[i].pos == sites[i - 1].pos && sites[i - 1].type == 2 )
		{
			if ( sites[i].type == 1 )	
				++sites[i].pos ;
		}
	}*/
}

void GradientDescentGammaDistribution( double &k, double &theta, double initK, double initTheta, double *x, double *z, int n ) 
{
	int i ;
	k = initK ;
	theta = initTheta ;
	double c = 0.00001 ;
	int iterCnt = 1 ;

	double sumZ = 0 ;
	double sumZX = 0 ;
	double sumZLogX = 0 ;
	double Hessian[2][2] ; // 0 for k, 1 for theta
	double inverseHessian[2][2] ;
	int tmp ;

	for ( i = 0 ; i < n ; ++i )	
	{
		sumZ += z[i] ;
		sumZX += z[i] * x[i]  ;
		sumZLogX += z[i] * log( x[i] ) ;
	}

	while ( 1 )
	{
		double gradK = 0 ;
		double gradTheta = 0 ;

		double prevK = k ;
		double prevTheta = theta ;

		double digammaK = digammal( k ) ;

		gradK = sumZ * ( -log( theta ) - digammaK ) + sumZLogX ;
		gradTheta = -sumZ * ( k / theta ) + sumZX / ( theta * theta ) ;

		Hessian[0][0] = -sumZ * trigamma( k, &tmp ) ;
		Hessian[0][1] = -sumZ / theta ; // \partial l / ( \partial k \partial theta)
		Hessian[1][0] = -sumZ / theta ;
		Hessian[1][1] = sumZ * k / ( theta * theta ) - 2 * sumZX / ( theta * theta * theta ) ;

		double det = Hessian[0][0] * Hessian[1][1] - Hessian[0][1] * Hessian[1][0] ;
		/*printf( "%lf %lf %lf %lf\n", sumZ, k, theta, sumZX ) ;	
		printf( "%lf %lf %lf\n", gradK, gradTheta, det ) ;	
		printf( "%lf %lf %lf %lf\n", Hessian[0][0], Hessian[0][1], Hessian[1][0], Hessian[1][1] ) ;*/
	
		if ( det == 0 )
		{
			k = k + c / iterCnt * gradK ;
			theta = theta + c / iterCnt * gradTheta ;
		}
		else
		{
			inverseHessian[0][0] = Hessian[1][1] / det ;
			inverseHessian[0][1] = -Hessian[0][1] / det ;
			inverseHessian[1][0] = -Hessian[1][0] / det ;
			inverseHessian[1][1] = Hessian[0][0] / det ;
			//printf( "%lf %lf %lf %lf: %lf\n=====\n", inverseHessian[0][0], inverseHessian[0][1], inverseHessian[1][0], inverseHessian[1][1],
			//	Hessian[1][0] * inverseHessian[0][1] + Hessian[1][1] * inverseHessian[1][1] ) ;	
			double step = 0.5 ;
			k = k - step * ( inverseHessian[0][0] * gradK + inverseHessian[0][1] * gradTheta ) ;
			theta = theta - step * ( inverseHessian[1][0] * gradK + inverseHessian[1][1] * gradTheta ) ;
		}

		if ( k <= 1e-6 )
			k = 1e-6 ;
		if ( theta <= 1e-6 )
			theta = 1e-6 ;

		double diff = ABS( prevK - k ) + ABS( prevTheta - theta ) ;
		if ( diff < 1e-5 ) 
			break ;
	
		++iterCnt ;
		if ( iterCnt == 1000 )
			break ;
	}
}


double MixtureGammaEM( double *x, int n, double &pi, double *k, double *theta, int iter = -1 )
{
	int i ;
	double *z = new double[n] ; // the expectation that it assigned to model 0.
	double *oneMinusZ = new double[n] ;

	int t = 0 ;
	while ( 1 )
	{
		double npi, nk[2], ntheta[2] ;
		double sum = 0 ;
		for ( i = 0 ; i < n ; ++i )	
		{
			//double lf0 = -k[0] * log( theta[0] ) + ( k[0] - 1 ) * log( cov[i]) - cov[i] / theta[0] - lgamma( k[0] ); 
			//double lf1 = -k[1] * log( theta[1] ) + ( k[1] - 1 ) * log( cov[i]) - cov[i] / theta[1] - lgamma( k[1] ); 
			//z[i] = exp( lf0 + log( pi ) ) / ( exp( lf0 + log( pi ) ) + exp( lf1 + log( 1 - pi ) ) ) ;
			z[i] = MixtureGammaAssignment( x[i], pi, k, theta ) ;
			oneMinusZ[i] = 1 - z[i] ;
			sum += z[i] ;
		}

		// compute new pi.
		npi = sum / n ;

		// Use gradient descent to compute new k and theta.
		GradientDescentGammaDistribution( nk[0], ntheta[0], k[0], theta[0], x, z, n ) ;
		GradientDescentGammaDistribution( nk[1], ntheta[1], k[1], theta[1], x, oneMinusZ,  n ) ;

		double diff ;
		diff = ABS( npi - pi ) + ABS( nk[0] - k[0] ) + ABS( nk[1] - k[1] )
			+ ABS( ntheta[0] - ntheta[0] ) + ABS( ntheta[1] - ntheta[1] ) ;
		if ( diff < 1e-4 )
			break ;
		pi = npi ;
		k[0] = nk[0] ;
		k[1] = nk[1] ;
		theta[0] = ntheta[0] ;
		theta[1] = ntheta[1] ;

		//printf( "%lf %lf %lf %lf %lf\n", pi, k[0], theta[0], k[1], theta[1] ) ;
		++t ;
		if ( iter != -1 && t >= iter )
			break ;
	}
	delete[] z ;
	delete[] oneMinusZ ;
	return 0 ;
}

int main( int argc, char *argv[] )
{
	int i, j ;

	bool noStats = false ;
	if ( argc < 3 )
	{
		fprintf( stderr, usage ) ;
		exit( 1 ) ;
	}

	gMinDepth = 2 ;

	for ( i = 3 ; i < argc ; ++i )
	{
		if ( !strcmp( argv[i], "--noStats" ) )
		{
			noStats = true ;
			continue ;
		}
		else if ( !strcmp( argv[i], "--minDepth" ) )
		{
			gMinDepth = atoi( argv[i + 1] ) ;
			++i ;
			continue ;
		}
		else
		{
			fprintf( stderr, "Unknown argument: %s\n", argv[i] ) ;
			return 0 ;
		}
	}

	Alignments alignments ;
	alignments.Open( argv[1] ) ;
	std::vector<struct _splitSite> splitSites ; // only compromised the 
	std::vector<struct _splitSite> allSplitSites ;

	// read in the splice site
	FILE *fp ;
	fp = fopen( argv[2], "r" ) ;
	char chrom[50] ;
	int64_t start, end ;
	int support ;
	char strand[3] ;
	int uniqSupport, secondarySupport, uniqEditDistance, secondaryEditDistance ;
	while ( fscanf( fp, "%s %" PRId64 " %" PRId64 " %d %s %d %d %d %d", chrom, &start, &end, &support, strand, 
				&uniqSupport, &secondarySupport, &uniqEditDistance, &secondaryEditDistance ) != EOF )
	{
		if ( support <= 0 )
			continue ;
		int chrId = alignments.GetChromIdFromName( chrom ) ; 
		struct _splitSite ss ;
		--start ;
		--end ;
		ss.pos = start ;
		ss.chrId = chrId ;
		ss.type = 2 ;
		ss.oppositePos = end ;
		splitSites.push_back( ss ) ;

		ss.pos = end ; 
		ss.type = 1 ;
		ss.oppositePos = start ;
		splitSites.push_back( ss ) ;
	}
	fclose( fp ) ;
	//printf( "ss:%d\n", splitSites.size() ) ;
	
	allSplitSites = splitSites ;

	CleanAndSortSplitSites( splitSites ) ;
	std::sort( allSplitSites.begin(), allSplitSites.end(), CompSplitSite ) ;
	//printf( "ss:%d\n", splitSites.size() ) ;

	// Build the blocks
	Blocks regions ;
	regions.BuildExonBlocks( alignments ) ;
	//printf( "%d\n", regions.exonBlocks.size() ) ;
	// Split the blocks using split site
	regions.SplitBlocks( alignments, splitSites ) ;
	//printf( "%d\n", regions.exonBlocks.size() ) ;

	// Recompute the coverage for each block. 
	alignments.Rewind() ;
	regions.ComputeDepth( alignments ) ;
	//printf( "%d\n", regions.exonBlocks.size() ) ;

	// Merge blocks that may have a hollow coverage by accident.
	regions.MergeNearBlocks() ;
	//printf( "%d\n", regions.exonBlocks.size() ) ;
	
	// Put the intron informations
	regions.AddIntronInformation( allSplitSites ) ;

	// Compute the average ratio against the left and right connected subexons.
	regions.ComputeRatios() ;
	
	if ( noStats ) 
	{ 
		// just output the subexons.
		if ( realpath( argv[1], buffer ) == NULL )
		{
			strcpy( buffer, argv[1] ) ;
		}
		printf( "#%s\n", buffer ) ;
		printf( "#fitted_ir_parameter_ratio: pi: -1 k0: -1 theta0: -1 k1: -1 theta1: -1\n" ) ;
		printf( "#fitted_ir_parameter_cov: pi: -1 k0: -1 theta0: -1 k1: -1 theta1: -1\n" ) ;
		
		int blockCnt = regions.exonBlocks.size() ;
		for ( int i = 0 ; i < blockCnt ; ++i )
		{
			struct _block &e = regions.exonBlocks[i] ;
			double avgDepth = (double)e.depthSum / ( e.end - e.start + 1 ) ;
			printf( "%s %" PRId64 " %" PRId64 " %d %d %lf -1 -1\n", alignments.GetChromName( e.chrId ), e.start + 1, e.end + 1, e.leftType, e.rightType, avgDepth ) ;
		}
		return 0 ;
	}

	// Extract the blocks for different events.
	int blockCnt = regions.exonBlocks.size() ;
	std::vector<struct _block> irBlocks ; // The regions corresponds to intron retention events.
	double *leftClassifier = new double[ blockCnt ] ; 
	double *rightClassifier = new double[ blockCnt ] ;
	std::vector<struct _block> overhangBlocks ; //blocks like (...[ or ]...)
	std::vector<struct _block> islandBlocks ; // blocks like (....) 

	for ( i = 0 ; i < blockCnt ; ++i )	
	{
		int idx = i ;
		
		int ltype = regions.exonBlocks[i].leftType ;
		int rtype = regions.exonBlocks[i].rightType ;
		leftClassifier[i] = -1 ;
		rightClassifier[i] = -1 ;
		
		double avgDepth = (double)regions.exonBlocks[i].depthSum / ( regions.exonBlocks[i].end - regions.exonBlocks[i].start + 1 ) ;
		
		if ( ltype == 2 && rtype == 1 )
		{
			// candidate intron retention.
			double ratio = regions.PickLeftAndRightRatio( regions.exonBlocks[i] ) ;

			//printf( "%lf %lf\n", regions.exonBlocks[i].leftRatio, regions.exonBlocks[i].rightRatio ) ;
			if ( ratio >= 0 )
			{
				regions.exonBlocks[i].ratio = ratio ;
				irBlocks.push_back( regions.exonBlocks[i] ) ;
				irBlocks[ irBlocks.size() - 1 ].contigId = i ;
			}
		}
		else if ( ( ltype == 0 && rtype == 1 ) || ( ltype == 2 && rtype == 0 ) )
		{
			// subexons like (...[ or ]...)
			double anchorAvg = 0 ;
			if ( ltype == 0 )
			{
				if ( idx < blockCnt - 1 && regions.exonBlocks[idx + 1].chrId == regions.exonBlocks[idx].chrId )
				{
					anchorAvg = regions.exonBlocks[idx+1].depthSum / ( regions.exonBlocks[idx+1].end - regions.exonBlocks[idx+1].start + 1 ) ;
				}
			}
			else if ( ltype == 2 )
			{
				
				if ( idx >= 1 && regions.exonBlocks[idx-1].chrId == regions.exonBlocks[idx].chrId )
				{
					anchorAvg = regions.exonBlocks[idx-1].depthSum / ( regions.exonBlocks[idx-1].end - regions.exonBlocks[idx-1].start + 1 ) ;
				}
			}
			if ( anchorAvg > 1 && avgDepth > 1 )
			{
				regions.exonBlocks[i].ratio = ( avgDepth - 1 ) / ( anchorAvg - 1 ) ;
				overhangBlocks.push_back( regions.exonBlocks[i] ) ;
				overhangBlocks[ overhangBlocks.size() - 1].contigId = i ;
			}
		}
		else if ( ltype == 0 && rtype == 0 )
		{
			islandBlocks.push_back( regions.exonBlocks[i] ) ;
			islandBlocks[ islandBlocks.size() - 1].contigId = i ;
		}
	}

	// Compute the histogram for each intron.
	int irBlockCnt = irBlocks.size() ;
	double *cov = new double[irBlockCnt] ;
	double *covRatio = new double[ irBlockCnt ] ;
	for ( i = 0 ; i < irBlockCnt ; ++i )
	{
		double avgDepth = regions.GetAvgDepth( irBlocks[i] ) ;
		//cov[i] = ( avgDepth - 1 ) / ( flankingAvg - 1 ) ;
		cov[i] = avgDepth - 1 ;

		covRatio[i] = irBlocks[i].leftRatio < irBlocks[i].rightRatio ? irBlocks[i].leftRatio : irBlocks[i].rightRatio ; 
		//cov[i] = avgDepth / anchorAvg ;
		//printf( "%"PRId64" %d %d: %lf %lf\n", irBlocks[i].depthSum, irBlocks[i].start, irBlocks[i].end, avgDepth, cov[i] ) ;
	}

	/*fp = fopen( "ratio.out", "r" ) ;
	int irBlockCnt = 0 ;
	double *cov = new double[10000] ;
	while ( 1 ) 
	{
		double r ;
		if ( fscanf( fp, "%lf", &r ) == EOF )
			break ;
		cov[ irBlockCnt ] = r ;
		++irBlockCnt ;
	}
	fclose( fp ) ;*/
	//printf( "hi\n" ) ;
	double piRatio = 0.8 ; // mixture coefficient for model 0 and 1
	double kRatio[2] = {1, 5} ;
	double thetaRatio[2] = {0.05, 0.1 } ;

	MixtureGammaEM( covRatio, irBlockCnt, piRatio, kRatio, thetaRatio ) ;

	// Test whether we should use single gamma distribution model or mixture model.
	// TODO: use a subset of all candidate to speed up this procedure.
	// Try fit a single gamma distribution model
	double sk, stheta ;
	double *all1 = new double[irBlockCnt] ;
	for ( i = 0 ; i < irBlockCnt ; ++i )
		all1[i] = 1 ;

	if ( piRatio >= 0.5 )
		GradientDescentGammaDistribution( sk, stheta, kRatio[0], thetaRatio[0], covRatio, all1, irBlockCnt ) ;
	else
		GradientDescentGammaDistribution( sk, stheta, kRatio[1], thetaRatio[1], covRatio, all1, irBlockCnt ) ;
	//sk = k[0] ;
	//stheta = theta[0] ;
	delete[] all1 ;
	//printf( "%lf %lf\n", sk, stheta ) ;
	double loglikelihoodMixture = 0 ;
	double loglikelihoodSingle = 0 ;
	for ( i = 0 ; i < irBlockCnt ; ++i )
	{
		double lf0 = LogGammaDensity( covRatio[i], kRatio[0], thetaRatio[0] ) ;
		double lf1 = LogGammaDensity( covRatio[i], kRatio[1], thetaRatio[1] ) ;
		double lfs = LogGammaDensity( covRatio[i], sk, stheta ) ;
		
		loglikelihoodMixture += log( exp( log( piRatio ) + lf0 ) + exp( log( 1 - piRatio ) + lf1 ) ) ;
		loglikelihoodSingle += lfs ;
	}
	//printf( "llmixture: %lf llsingle: %lf\n", loglikelihoodMixture, loglikelihoodSingle ) ;
	/*int simulateTimes = 1000 ;
	double *simulateLogDiff = new double[simulateTimes] ;
	std::mt19937 generator(1701) ;
	std::gamma_distribution<double> gammaDistn( sk, stheta ) ;

	for ( j = 0 ; j < simulateTimes ; ++j )
	{
		double llmixture = 0 ;
		double llsingle = 0 ;
		for ( i = 0 ; i < irBlockCnt ; ++i )			
		{
			double x = gammaDistn( generator ) ;
			llmixture += log( exp( log( pi ) + LogGammaDensity( x, k[0], theta[0] ) ) + exp( log( 1 - pi ) + LogGammaDensity( x, k[1], theta[1] ) ) ) ;
			llsingle += LogGammaDensity( x, sk, stheta ) ;
		}
		printf( "%lf %lf\n", llmixture, llsingle ); 
		simulateLogDiff[j] = llmixture - llsingle ;
	}
	qsort( simulateLogDiff, simulateTimes, sizeof( double ), CompDouble ) ;

	double testDiff = loglikelihoodMixture - loglikelihoodSingle ;
	printf( "%lf %lf\n", simulateLogDiff[int( 0.95 * simulateTimes )], testDiff ) ;
	delete []simulateLogDiff ;*/

	double bicMixture = 2 * log( irBlockCnt ) * 5 - 2 * loglikelihoodMixture ;
	double bicSingle = 2 * log( irBlockCnt ) * 2 - 2 * loglikelihoodSingle ;
	//printf( "bic: %lf %lf\n", bicMixture, bicSingle ) ;
	
	//double aicMixture = 2 * 5 - 2 * loglikelihoodMixture ;
	//double aicSingle = 2 * 2 - 2 * loglikelihoodSingle ;
	//printf( "aic: %lf %lf\n", aicMixture, aicSingle ) ;
	
	double piCov = piRatio ; // mixture coefficient for model 0 and 1
	double kCov[2] = {1, 5} ;
	double thetaCov[2] = {3, 3 } ;
	
	if ( bicSingle < bicMixture )
	{
		piRatio = 1 ;
		kRatio[0] = sk ;
		thetaRatio[0] = stheta ;
		kRatio[1] = 0 ;
		thetaRatio[1] = 0 ;

		piCov =1 ;
		kCov[0] = 0 ;
		thetaCov[0] = 0 ;
		kCov[1] = 0 ;
		thetaCov[1] = 0 ;
	}
	else
	{
		// only do one iteration of EM, so that pi does not change.
		MixtureGammaEM( cov, irBlockCnt, piCov, kCov, thetaCov, 1 ) ;	
		piCov = piRatio ;
	}

	for ( i = 0 ; i < irBlockCnt ; ++i )
	{
		//double lf0 = -k[0] * log( theta[0] ) + ( k[0] - 1 ) * log( cov[i]) - cov[i] / theta[0] - lgamma( k[0] ) ;
		//double lf1 = -k[1] * log( theta[1] ) + ( k[1] - 1 ) * log( cov[i]) - cov[i] / theta[1] - lgamma( k[1] ) ;

		double p1, p2, p ;
		if ( piRatio != 1 )
		{
			p1 = MixtureGammaAssignment( covRatio[i], piRatio, kRatio, thetaRatio ) ;
			p2 = MixtureGammaAssignment( cov[i], piCov, kCov, thetaCov  ) ;
			p = p1>p2 ? p1 : p2 ;
		}
		else
			p = -1 ;
		//printf( "%d %d: avg: %lf ratio: %lf p: %lf\n", irBlocks[i].start, irBlocks[i].end, irBlocks[i].depthSum / (double)( irBlocks[i].end - irBlocks[i].start + 1 ), covRatio[i], 
		//		p ) ;	
		leftClassifier[ irBlocks[i].contigId ] = p ;
		rightClassifier[ irBlocks[i].contigId ] = p ;
	}

	// Process the classifier for overhang subexons and the subexons to see whether we need soft boundary
	for ( i = 0 ; i < blockCnt ; ++i )		
	{
		struct _block &e = regions.exonBlocks[i] ;
		if ( ( e.leftType == 0 && e.rightType == 1 ) || 
			e.leftType == 1 )
		{
			if ( 2 * e.leftRatio >= 0 )
				leftClassifier[i] = 2 * alnorm( 2 * e.leftRatio, true ) ;
			else
				leftClassifier[i] = 1 ;
		}
		if ( ( e.rightType == 0 && e.leftType == 2 ) || 
			e.rightType == 2 )
		{
			if ( 2 * e.rightRatio >= 0 )
				rightClassifier[i] = 2 * alnorm( 2 * e.rightRatio, true ) ;
			else
				rightClassifier[i] = 1 ;
		}
	}

	// Process the result for subexons seems like single-exon transcript (...)
	int islandBlockCnt = islandBlocks.size() ;
	for ( i = 0, j = 0 ; i < islandBlockCnt ; ++i )
	{
		if ( regions.GetAvgDepth( islandBlocks[i] ) != regions.GetAvgDepth( islandBlocks[j] ) )
			++j ;
		leftClassifier[ islandBlocks[i].contigId ] = 1 - (j + 1) / (double)( islandBlockCnt + 1 ) ;
		rightClassifier[ islandBlocks[i].contigId ] = 1 - (j + 1) / (double)( islandBlockCnt + 1 ) ;
	}

	
	// Output the result.
	if ( realpath( argv[1], buffer ) == NULL )
	{
		strcpy( buffer, argv[1] ) ;
	}
	printf( "#%s\n", buffer ) ;
	printf( "#fitted_ir_parameter_ratio: pi: %lf k0: %lf theta0: %lf k1: %lf theta1: %lf\n", piRatio, kRatio[0], thetaRatio[0], kRatio[1], thetaRatio[1] ) ;
	printf( "#fitted_ir_parameter_cov: pi: %lf k0: %lf theta0: %lf k1: %lf theta1: %lf\n", piCov, kCov[0], thetaCov[0], kCov[1], thetaCov[1] ) ;

	for ( int i = 0 ; i < blockCnt ; ++i )
	{
		struct _block &e = regions.exonBlocks[i] ;
		double avgDepth = regions.GetAvgDepth( e ) ;
		printf( "%s %" PRId64 " %" PRId64 " %d %d %lf %lf %lf %lf %lf ", alignments.GetChromName( e.chrId ), e.start + 1, e.end + 1, e.leftType, e.rightType, avgDepth, 
			e.leftRatio, e.rightRatio, leftClassifier[i], rightClassifier[i] ) ;
		int prevCnt = e.prevCnt ;
		if ( i > 0 && e.start == regions.exonBlocks[i - 1].end + 1 &&
			e.leftType == regions.exonBlocks[i - 1].rightType )
		{
			printf( "%d ", prevCnt + 1 ) ;
			for ( j = 0 ; j < prevCnt ; ++j )
				printf( "%" PRId64 " ", regions.exonBlocks[ e.prev[j] ].end + 1 ) ;
			printf( "%" PRId64 " ", regions.exonBlocks[i - 1].end + 1 ) ;
		}
		else
		{
			printf( "%d ", prevCnt ) ;
			for ( j = 0 ; j < prevCnt ; ++j )
				printf( "%" PRId64 " ", regions.exonBlocks[ e.prev[j] ].end + 1 ) ;
		}

		int nextCnt = e.nextCnt ;
		if ( i < blockCnt - 1 && e.end == regions.exonBlocks[i + 1].start - 1 &&
			e.rightType == regions.exonBlocks[i + 1].leftType )
		{
			printf( "%d %" PRId64 " ", nextCnt + 1, regions.exonBlocks[i + 1].start + 1 ) ;
		}
		else
			printf( "%d ", nextCnt ) ;
		for ( j = 0 ; j < nextCnt ; ++j )
			printf( "%" PRId64 " ", regions.exonBlocks[ e.next[j] ].start + 1 ) ;
		printf( "\n" ) ;
	}

	delete[] cov ;
	delete[] covRatio ;
	delete[] leftClassifier ;
	delete[] rightClassifier ;
}

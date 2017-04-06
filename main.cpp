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

#include <random>

#include "alignments.hpp"
#include "blocks.hpp"
#include "gamma.hpp"

#define ABS(x) ((x)<0?-(x):(x))

char usage[] = "./subexon-info alignment.bam intron.splice [options]\n" ;
char buffer[4096] ;

bool CompSplitSite( struct _splitSite a, struct _splitSite b )
{
	if ( a.chrId < b.chrId )
		return true ;
	else if ( a.chrId > b.chrId )
		return false ;
	else if ( a.pos != b.pos )
		return a.pos < b.pos ;
	else
		return a.type < b.type ; // We want the start of exons comes first, since we are scanning the genome from left to right,
					 // so we can terminate the extension early and create a single-base subexon later.
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

	if ( argc != 3 )
	{
		fprintf( stderr, "./a.out alignments.bam splice_site\n" ) ;
		exit( 1 ) ;
	}

	Alignments alignments ;
	alignments.Open( argv[1] ) ;
	std::vector<struct _splitSite> splitSites ;

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
		splitSites.push_back( ss ) ;

		ss.pos = end ; 
		ss.type = 1 ;
		splitSites.push_back( ss ) ;
	}
	fclose( fp ) ;
	//printf( "ss:%d\n", splitSites.size() ) ;
	CleanAndSortSplitSites( splitSites ) ;
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

	// Compute the possible intron-retentioned blocks.
	int blockCnt = regions.exonBlocks.size() ;
	std::vector<struct _block> irBlocks ; // The regions corresponds to intron retention events.
	std::vector<int> irBlockIdx ; // their index in the exonBlocks list
	double *irClassifier = new double[ blockCnt ] ;
	for ( i = 0 ; i < blockCnt ; ++i )	
	{
		double flankingAvg = 0 ;
		double flankingLen = 0 ;
		int idx = i ;
		double anchorAvg = 0 ;
		if ( idx >= 1 && regions.exonBlocks[idx-1].chrId == regions.exonBlocks[idx].chrId )
		{
			flankingAvg += regions.exonBlocks[idx-1].depthSum ;
			flankingLen += ( regions.exonBlocks[idx-1].end - regions.exonBlocks[idx-1].start + 1 ) ;
			//printf( "left %d: %d %d %"PRId64"\n", idx - 1, regions.exonBlocks[idx-1].start, regions.exonBlocks[idx-1].end, regions.exonBlocks[idx-1].depthSum) ;
			anchorAvg = regions.exonBlocks[idx-1].depthSum / ( regions.exonBlocks[idx-1].end - regions.exonBlocks[idx-1].start + 1 ) ;
		}
		if ( idx < blockCnt - 1 && regions.exonBlocks[idx + 1].chrId == regions.exonBlocks[idx].chrId )
		{
			flankingAvg += regions.exonBlocks[idx+1].depthSum ;
			flankingLen += ( regions.exonBlocks[idx+1].end - regions.exonBlocks[idx+1].start + 1 ) ;
			//printf( "right %d: %d %d %"PRId64"\n", idx + 1, regions.exonBlocks[idx+1].start, regions.exonBlocks[idx+1].end, regions.exonBlocks[idx+1].depthSum) ;
			double tmp = regions.exonBlocks[idx+1].depthSum / ( regions.exonBlocks[idx+1].end - regions.exonBlocks[idx+1].start + 1 ) ;
			if ( tmp > anchorAvg )
				anchorAvg = tmp ;
		}
		flankingAvg /= flankingLen ;
		double avgDepth = (double)regions.exonBlocks[i].depthSum / ( regions.exonBlocks[i].end - regions.exonBlocks[i].start + 1 ) ;
		
		irClassifier[i] = -1 ;
		if ( anchorAvg > 1 && avgDepth > 1 && regions.exonBlocks[i].leftType == 2 && regions.exonBlocks[i].rightType == 1 )
		{
			irBlocks.push_back( regions.exonBlocks[i] ) ;
			irBlockIdx.push_back( i ) ;
		}
	}
	//printf( "%d\n", irBlocks.size() ) ;
	// Compute the histogram for each intron.
	int irBlockCnt = irBlocks.size() ;
	double *cov = new double[irBlockCnt] ;
	double *covRatio = new double[ irBlockCnt ] ;
	for ( i = 0 ; i < irBlockCnt ; ++i )
	{
		double flankingAvg = 0 ;
		double flankingLen = 0 ;
		int idx = irBlockIdx[i] ;

		double anchorAvg = 0 ;	
		double avgDepth = (double)irBlocks[i].depthSum / (irBlocks[i].end - irBlocks[i].start + 1 ) ;
		cov[i] = 10000000 ;
		if ( idx >= 1 && regions.exonBlocks[idx-1].chrId == regions.exonBlocks[idx].chrId )
		{
			flankingAvg += regions.exonBlocks[idx-1].depthSum ;
			flankingLen += ( regions.exonBlocks[idx-1].end - regions.exonBlocks[idx-1].start + 1 ) ;
			//printf( "left %d: %d %d %"PRId64"\n", idx - 1, regions.exonBlocks[idx-1].start, regions.exonBlocks[idx-1].end, regions.exonBlocks[idx-1].depthSum) ;
			anchorAvg = regions.exonBlocks[idx-1].depthSum / ( regions.exonBlocks[idx-1].end - regions.exonBlocks[idx-1].start + 1 ) ;
		}
		if ( idx < blockCnt - 1 && regions.exonBlocks[idx + 1].chrId == regions.exonBlocks[idx].chrId )
		{
			flankingAvg += regions.exonBlocks[idx+1].depthSum ;
			flankingLen += ( regions.exonBlocks[idx+1].end - regions.exonBlocks[idx+1].start + 1 ) ;
			//printf( "right %d: %d %d %"PRId64"\n", idx + 1, regions.exonBlocks[idx+1].start, regions.exonBlocks[idx+1].end, regions.exonBlocks[idx+1].depthSum) ;
			
			double tmp = regions.exonBlocks[idx+1].depthSum / ( regions.exonBlocks[idx+1].end - regions.exonBlocks[idx+1].start + 1 ) ;
			if ( tmp > anchorAvg )
				anchorAvg = tmp ;
		}
		flankingAvg /= flankingLen ;
		//cov[i] = ( avgDepth - 1 ) / ( flankingAvg - 1 ) ;
		cov[i] = avgDepth - 1 ;
		covRatio[i] = ( avgDepth - 1 ) / ( anchorAvg - 1 ) ;
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


	// Output the result.
	if ( realpath( argv[1], buffer ) == NULL )
	{
		strcpy( buffer, argv[1] ) ;
	}
	printf( "#%s\n", buffer ) ;
	printf( "#fitted_ir_parameter_ratio: pi: %lf k0: %lf theta0: %lf k1: %lf theta1: %lf\n", piRatio, kRatio[0], thetaRatio[0], kRatio[1], thetaRatio[1] ) ;
	printf( "#fitted_ir_parameter_cov: pi: %lf k0: %lf theta0: %lf k1: %lf theta1: %lf\n", piCov, kCov[0], thetaCov[0], kCov[1], thetaCov[1] ) ;
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
		irClassifier[ irBlockIdx[i] ] = p ;
	}

	for ( int i = 0 ; i < blockCnt ; ++i )
	{
		struct _block &e = regions.exonBlocks[i] ;
		double avgDepth = (double)e.depthSum / ( e.end - e.start + 1 ) ;
		printf( "%s %" PRId64 " %" PRId64 " %d %d %lf %lf\n", alignments.GetChromName( e.chrId ), e.start + 1, e.end + 1, e.leftType, e.rightType, avgDepth, irClassifier[i] ) ;
	}

	delete[] cov ;
	delete[] covRatio ;
	delete[] irClassifier ;
}

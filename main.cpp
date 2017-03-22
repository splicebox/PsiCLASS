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
#include "gamma.hpp"

#define ABS(x) ((x)<0?-(x):(x))

bool CompSplitSite( struct _splitSite a, struct _splitSite b )
{
	if ( a.chrId < b.chrId )
		return true ;
	else if ( a.chrId > b.chrId )
		return false ;
	else if ( a.pos != b.pos )
		return a.pos < b.pos ;
	else
		return b.type - a.type ; // We want the end of exons comes first, since we are scanning the genome from left to right.
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

double MixtureGammaAssignment( double x, double pi, double* k, double *theta )
{
	double lf0 = -k[0] * log( theta[0] ) + ( k[0] - 1 ) * log( x ) - x / theta[0] - lgamma( k[0] ) ;
	double lf1 = -k[1] * log( theta[1] ) + ( k[1] - 1 ) * log( x ) - x / theta[1] - lgamma( k[1] ) ;

	return (double)1.0 / ( 1.0 + exp( lf1 + log( 1 - pi ) - lf0 - log( pi ) ) ) ;

}


int main( int argc, char *argv[] )
{
	int i ;

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
	while ( fscanf( fp, "%s %"PRId64" %"PRId64" %d %s %d %d %d %d", chrom, &start, &end, &support, strand, 
				&uniqSupport, &secondarySupport, &uniqEditDistance, &secondaryEditDistance ) != EOF )
	{
		int chrId = alignments.GetChromIdFromName( chrom ) ; 	
		struct _splitSite ss ;
		--start ;
		--end ;
		ss.pos = start ;
		ss.chrId = chrId ;
		ss.type = 2 ;
		splitSites.push_back( ss ) ;

		ss.pos = end - 1 ; // Note the adjust to the right-end.
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
	for ( i = 0 ; i < blockCnt ; ++i )	
	{
		double flankingAvg = 0 ;
		double flankingLen = 0 ;
		int idx = i ;
		double anchorAvg = 0 ;
		if ( idx >= 1 )
		{
			flankingAvg += regions.exonBlocks[idx-1].depthSum ;
			flankingLen += ( regions.exonBlocks[idx-1].end - regions.exonBlocks[idx-1].start + 1 ) ;
			//printf( "left %d: %d %d %"PRId64"\n", idx - 1, regions.exonBlocks[idx-1].start, regions.exonBlocks[idx-1].end, regions.exonBlocks[idx-1].depthSum) ;
			anchorAvg = regions.exonBlocks[idx-1].depthSum / ( regions.exonBlocks[idx-1].end - regions.exonBlocks[idx-1].start + 1 ) ;
		}
		if ( idx < blockCnt - 1 )
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
	for ( i = 0 ; i < irBlockCnt ; ++i )
	{
		double flankingAvg = 0 ;
		double flankingLen = 0 ;
		int idx = irBlockIdx[i] ;

		double anchorAvg = 0 ;	
		double avgDepth = (double)irBlocks[i].depthSum / (irBlocks[i].end - irBlocks[i].start + 1 ) ;
		cov[i] = 10000000 ;
		if ( idx >= 1 )
		{
			flankingAvg += regions.exonBlocks[idx-1].depthSum ;
			flankingLen += ( regions.exonBlocks[idx-1].end - regions.exonBlocks[idx-1].start + 1 ) ;
			//printf( "left %d: %d %d %"PRId64"\n", idx - 1, regions.exonBlocks[idx-1].start, regions.exonBlocks[idx-1].end, regions.exonBlocks[idx-1].depthSum) ;
			anchorAvg = regions.exonBlocks[idx-1].depthSum / ( regions.exonBlocks[idx-1].end - regions.exonBlocks[idx-1].start + 1 ) ;
		}
		if ( idx < blockCnt - 1 )
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
		cov[i] = ( avgDepth - 1 ) / ( anchorAvg - 1 ) ;
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
	double pi = 0.8 ; // mixture coefficient for model 0 and 1
	double k[2] = {1, 5} ;
	double theta[2] = {0.05, 0.1 } ;

	double *z = new double[irBlockCnt] ; // the expectation that it assigned to model 0.
	double *oneMinusZ = new double[irBlockCnt] ;
	while ( 1 )
	{
		double npi, nk[2], ntheta[2] ;
		double sum = 0 ;
		for ( i = 0 ; i < irBlockCnt ; ++i )	
		{
			//double lf0 = -k[0] * log( theta[0] ) + ( k[0] - 1 ) * log( cov[i]) - cov[i] / theta[0] - lgamma( k[0] ); 
			//double lf1 = -k[1] * log( theta[1] ) + ( k[1] - 1 ) * log( cov[i]) - cov[i] / theta[1] - lgamma( k[1] ); 
			//z[i] = exp( lf0 + log( pi ) ) / ( exp( lf0 + log( pi ) ) + exp( lf1 + log( 1 - pi ) ) ) ;
			z[i] = MixtureGammaAssignment( cov[i], pi, k, theta ) ;
			oneMinusZ[i] = 1 - z[i] ;
			sum += z[i] ;
		}

		// compute new pi.
		npi = sum / irBlockCnt ;

		// Use gradient descent to compute new k and theta.
		GradientDescentGammaDistribution( nk[0], ntheta[0], k[0], theta[0], cov, z, irBlockCnt ) ;
		GradientDescentGammaDistribution( nk[1], ntheta[1], k[1], theta[1], cov, oneMinusZ,  irBlockCnt ) ;

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
	}
	printf( "fitted result: pi=%lf k0=%lf theta0=%lf k1=%lf theta1=%lf\n", pi, k[0], theta[0], k[1], theta[1] ) ;
	// Output
	for ( i = 0 ; i < irBlockCnt ; ++i )
	{
		//double lf0 = -k[0] * log( theta[0] ) + ( k[0] - 1 ) * log( cov[i]) - cov[i] / theta[0] - lgamma( k[0] ) ;
		//double lf1 = -k[1] * log( theta[1] ) + ( k[1] - 1 ) * log( cov[i]) - cov[i] / theta[1] - lgamma( k[1] ) ;
		
		printf( "%d %d: avg: %lf ratio: %lf p: %lf\n", irBlocks[i].start, irBlocks[i].end, irBlocks[i].depthSum / (double)( irBlocks[i].end - irBlocks[i].start + 1 ), cov[i], 
				MixtureGammaAssignment( cov[i], pi, k, theta ) ) ;		
	}

	/*for ( int i = 0 ; i < blockCnt ; ++i )
	{
		struct _block &e = regions.exonBlocks[i] ;
		//printf( "%s %"PRId64" %"PRId64" %"PRId64"\n", alignments.GetChromName( e.chrId ), e.start + 1, e.end + 1, e.depthSum ) ;
	}*/

	delete[] z ;
	delete[] cov ;
}

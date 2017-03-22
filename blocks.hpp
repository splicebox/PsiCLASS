// The class manage the blocks
// Li Song

#ifndef _LSONG_RSCAF_BLOCKS_HEADER
#define _LSONG_RSCAF_BLOCKS_HEADER

#include <stdlib.h> 
#include <vector>
#include <map>
#include <assert.h>
#include <math.h>
#include <set>
#include <inttypes.h>


#include "defs.h"

extern int minimumSupport ;
extern int minimumEffectiveLength ;
extern int kmerSize ;
extern bool VERBOSE ;
extern FILE *fpOut ;
extern bool aggressiveMode ;

struct _splitSite // means the next position belongs to another block
{
	int64_t pos ;
	//int strand ;
	int chrId ;
	int type ; // 1-start of an exon. 2-end of an exon.
} ;

struct _block
{
	int chrId ;
	int contigId ;
	int64_t start, end ;
	int64_t leftSplice, rightSplice ; // record the leftmost and rightmost coordinates of a splice site within this block.
	int64_t depthSum ;

	int leftType, rightType ; //0-soft boundary, 1-start of an exon, 2-end of an exon.
} ;

class Blocks
{
	private:
		std::map<int, int> exonBlocksChrIdOffset ;

		int64_t Overlap( int64_t s0, int64_t e0, int64_t s1, int64_t e1 )
		{
			if ( e0 < s1 || s0 > e1 )
				return 0 ;
			int64_t s, e ;
			s = s0 > s1 ? s0 : s1 ;
			e = e0 < e1 ? e0 : e1 ;
			return e - s + 1 ;
		}

		void Split( const char *s, char delimit, std::vector<std::string> &fields )
		{
			int i ;
			fields.clear() ;
			if ( s == NULL )
				return ;

			std::string f ;
			for ( i = 0 ; s[i] ; ++i )
			{
				if ( s[i] == delimit || s[i] == '\n' )	
				{
					fields.push_back( f ) ;
					f.clear() ;
				}
				else
					f.append( 1, s[i] ) ;
			}
			fields.push_back( f ) ;
			f.clear() ;
		}

	public:
		std::vector<struct _block> exonBlocks ;

		Blocks() { } 	
		~Blocks() 
		{
		}

		int BuildExonBlocks( Alignments &alignments )
		{
			unsigned int tag = 0 ;
			while ( alignments.Next() )
			{
				int i, j, k ;
				int segCnt = alignments.segCnt ;
				struct _pair *segments = alignments.segments ;

				while ( tag < exonBlocks.size() && ( exonBlocks[tag].end < segments[0].a - 1 
							|| exonBlocks[tag].chrId != alignments.GetChromId() ) )  
				{
					++tag ;
				}

				for ( i = 0 ; i < segCnt ; ++i )
				{
					//if ( i == 0 )
					//	printf( "hi %d %d %d\n", i, segments[i].a, segments[i].b ) ;
					for ( j = tag ; j < (int)exonBlocks.size() ; ++j )
					{
						if ( exonBlocks[j].end >= segments[i].a - 1 )
							break ;
					}

					if ( j >= (int)exonBlocks.size() )
					{
						// Append a new block
						struct _block newSeg ;
						newSeg.chrId = alignments.GetChromId() ;
						newSeg.start = segments[i].a ;
						newSeg.end = segments[i].b ;

						newSeg.leftSplice = -1 ;
						newSeg.rightSplice = -1 ;
						if ( i > 0 )
							newSeg.leftSplice = segments[i].a ; 
						if ( i < segCnt - 1 )
							newSeg.rightSplice = segments[i].b ;

						exonBlocks.push_back( newSeg ) ;
					}
					else if ( exonBlocks[j].end < segments[i].b || 
							( exonBlocks[j].start > segments[i].a && exonBlocks[j].start <= segments[i].b + 1 ) ) 
					{
						// If overlaps with a current exon block, so we extend it
						if ( exonBlocks[j].end < segments[i].b )
						{
							// extends toward right 
							exonBlocks[j].end = segments[i].b ;
							if ( i > 0 && ( exonBlocks[j].leftSplice == -1 || segments[i].a < exonBlocks[j].leftSplice ) )
								exonBlocks[j].leftSplice = segments[i].a ;
							if ( i < segCnt - 1 && segments[i].b > exonBlocks[j].rightSplice )
								exonBlocks[j].rightSplice = segments[i].b ;

							// Merge with next few exon blocks
							for ( k = j + 1 ; k < (int)exonBlocks.size() ; ++k )
							{
								if ( exonBlocks[k].start <= exonBlocks[j].end + 1 )
								{
									if ( exonBlocks[k].end > exonBlocks[j].end )
										exonBlocks[j].end = exonBlocks[k].end ;

									if ( exonBlocks[k].leftSplice != -1 && ( exonBlocks[j].leftSplice == -1 || exonBlocks[k].leftSplice < exonBlocks[j].leftSplice ) )
										exonBlocks[j].leftSplice = exonBlocks[k].leftSplice ;
									if ( exonBlocks[k].rightSplice != -1 && exonBlocks[k].rightSplice > exonBlocks[j].rightSplice )
										exonBlocks[j].rightSplice = exonBlocks[k].rightSplice ;

								}
								else
									break ;
							}

							if ( k > j + 1 )
							{
								// Remove the merged blocks
								int a, b ;
								for ( a = j + 1, b = k ; b < (int)exonBlocks.size() ; ++a, ++b )
									exonBlocks[a] = exonBlocks[b] ;
								for ( a = 0 ; a < k - ( j + 1 ) ; ++a )
									exonBlocks.pop_back() ;
							}
						}
						else if ( exonBlocks[j].start > segments[i].a && exonBlocks[j].start <= segments[i].b + 1 ) 
						{
							// extends toward left
							exonBlocks[j].start = segments[i].a ;
							if ( i > 0 && ( exonBlocks[j].leftSplice == -1 || segments[i].a < exonBlocks[j].leftSplice ) )
								exonBlocks[j].leftSplice = segments[i].a ;
							if ( i < segCnt - 1 && segments[i].b > exonBlocks[j].rightSplice )
								exonBlocks[j].rightSplice = segments[i].b ;

							// Merge with few previous exon blocks
							for ( k = j - 1 ; k >= 0 ; --k )
							{
								if ( exonBlocks[k].end >= exonBlocks[k + 1].start - 1 )
								{
									if ( exonBlocks[k + 1].start < exonBlocks[k].start )
									{
										exonBlocks[k].start = exonBlocks[k + 1].start ;
									}

									if ( exonBlocks[k].leftSplice != -1 && ( exonBlocks[j].leftSplice == -1 || exonBlocks[k].leftSplice < exonBlocks[j].leftSplice ) )
										exonBlocks[j].leftSplice = exonBlocks[k].leftSplice ;
									if ( exonBlocks[k].rightSplice != -1 && exonBlocks[k].rightSplice > exonBlocks[j].rightSplice )
										exonBlocks[j].rightSplice = exonBlocks[k].rightSplice ;

								}
								else
									break ;
							}

							if ( k < j - 1 )
							{
								int a, b ;
								for ( a = k + 2, b = j + 1 ; b < (int)exonBlocks.size() ; ++a, ++b )
									exonBlocks[a] = exonBlocks[b] ;
								for ( a = 0 ; a < ( j - 1 ) - k ; ++a )
									exonBlocks.pop_back() ;

							}
						}
					}
					else if ( exonBlocks[j].start > segments[i].b + 1 )
					{
						int size = exonBlocks.size() ;
						int a ;
						// No overlap, insert a new block
						struct _block newSeg ;
						newSeg.chrId = alignments.GetChromId() ;
						newSeg.start = segments[i].a ;
						newSeg.end = segments[i].b ;

						newSeg.leftSplice = -1 ;
						newSeg.rightSplice = -1 ;
						if ( i > 0 )
							newSeg.leftSplice = segments[i].a ; 
						if ( i < segCnt - 1 )
							newSeg.rightSplice = segments[i].b ;

						// Insert at position j
						exonBlocks.push_back( newSeg ) ;	
						for ( a = size ; a > j ; --a )
							exonBlocks[a] = exonBlocks[a - 1] ;
						exonBlocks[a] = newSeg ;
					}
					else
					{
						// The segment is contained in j
					}
				}
			}

			/*for ( int i = 0 ; i < (int)exonBlocks.size() ; ++i )
			  {
			  printf( "%d %d %d\n", exonBlocks[i].start, exonBlocks[i].end, exonBlocks[i].support.GetCount() ) ;
			  }*/	

			if ( exonBlocks.size() > 0 )
			{
				// Build the map for the offsets of chr id in the exonBlock list.
				exonBlocksChrIdOffset[ exonBlocks[0].chrId] = 0 ;
				int cnt = exonBlocks.size() ;
				for ( int i = 1 ; i < cnt ; ++i )
				{
					if ( exonBlocks[i].chrId != exonBlocks[i - 1].chrId )
						exonBlocksChrIdOffset[ exonBlocks[i].chrId ] = i ;
				}

				for ( int i = 0 ; i < cnt ; ++i )
				{
					exonBlocks[i].contigId = exonBlocks[i].chrId ;

					exonBlocks[i].leftType = 0 ;
					exonBlocks[i].rightType = 0 ;
				}
			}
			return exonBlocks.size() ;
		}


		void SplitBlocks( Alignments &alignments, std::vector< struct _splitSite > &splitSites )	
		{
			std::vector<struct _block> rawExonBlocks = exonBlocks ;
			int i, j ;
			int tag = 0 ;
			int bsize = rawExonBlocks.size() ; 
			int ssize = splitSites.size() ;

			// Make sure not overflow
			struct _splitSite tmp ;
			tmp.pos = -1 ;
			splitSites.push_back( tmp ) ;

			exonBlocks.clear() ;
			for ( i = 0 ; i < bsize ; ++i )
			{
				while ( tag < ssize && ( splitSites[tag].chrId < rawExonBlocks[i].chrId ||
							( splitSites[tag].chrId == rawExonBlocks[i].chrId && splitSites[tag].pos < rawExonBlocks[i].start ) ) )
					++tag ;

				int l ;
				for ( l = tag ; l < ssize ; ++l )
					if ( splitSites[l].chrId != rawExonBlocks[i].chrId || splitSites[l].pos > rawExonBlocks[i].end )
						break ;
				int64_t start = rawExonBlocks[i].start ;
				int64_t end = rawExonBlocks[i].end ;
				//printf( "%lld %lld\n", start, end ) ;
				for ( j = tag ; j <= l ; ++j )
				{
					int leftType = 0 ;
					int rightType = 0 ;
					if ( j > tag )	
					{
						start = end + 1 ; // end is from previous stage
						leftType = splitSites[j - 1].type ;
					}
					else
						start = rawExonBlocks[i].start ;

					if ( j <= l - 1 )
					{
						end = splitSites[j].pos ;
						rightType = splitSites[j].type ;
					}
					else
						end = rawExonBlocks[i].end ;

					struct _block tmpB ;
					tmpB = rawExonBlocks[i] ;
					tmpB.start = start ;
					tmpB.end = end ;
					tmpB.depthSum = 0 ;
					//printf( "\t%lld %lld\n", start, end ) ;
					tmpB.leftType = leftType ;
					tmpB.rightType = rightType ;
					exonBlocks.push_back( tmpB ) ;
					// handle the soft boundary is the same as the hard boundary case 
					// or adjacent splice sites
					/*if ( j == tag && start == end )
					{
						exonBlocks.pop_back() ;
						--end ;
					}
					else if ( j == l && start > end )
					{
						exonBlocks.pop_back() ;
					}*/
					if ( start >= end )
					{
						exonBlocks.pop_back() ;
						--end ;
					}
				}
			}
		}

		void ComputeDepth( Alignments &alignments ) 
		{
			// Go through the alignment list again to fill in the depthSum;
			int i ;
			int j ;
			int tag = 0 ;
			int blockCnt = exonBlocks.size() ;
			while ( alignments.Next() )
			{
				int segCnt = alignments.segCnt ;
				struct _pair *segments = alignments.segments ;

				while ( tag < blockCnt && ( exonBlocks[tag].chrId < alignments.GetChromId() || 
							( exonBlocks[tag].chrId == alignments.GetChromId() && exonBlocks[tag].end < segments[0].a ) ) )
					++tag ;
				for ( i = 0 ; i < segCnt ; ++i )
				{
					for ( j = tag ; j < blockCnt && ( exonBlocks[j].chrId == alignments.GetChromId() && exonBlocks[j].start <= segments[i].b ) ; ++j )
					{
						exonBlocks[j].depthSum += Overlap( segments[i].a, segments[i].b, exonBlocks[j].start, exonBlocks[j].end ) ; 
					}
				}
			}
		}

		// If two blocks whose soft boundary are close to each other, we can merge them.
		void MergeNearBlocks()
		{
			std::vector<struct _block> rawExonBlocks = exonBlocks ;
			int i, j, k ;
			int tag = 0 ;
			int bsize = rawExonBlocks.size() ; 
			
			exonBlocks.clear() ;
			exonBlocks.push_back( rawExonBlocks[0] ) ;
			k = 0 ;
			for ( i = 1 ; i < bsize ; ++i )	
			{
				if ( rawExonBlocks[i].chrId == exonBlocks[k].chrId 
					&& rawExonBlocks[i].leftType == 0 && exonBlocks[k].rightType == 0 
					&& rawExonBlocks[i].start - exonBlocks[k].end - 1 <= 15 )
				{
					exonBlocks[k].end = rawExonBlocks[i].end ;
					exonBlocks[k].rightType = rawExonBlocks[i].rightType ;
					exonBlocks[k].depthSum += rawExonBlocks[i].depthSum ;
					if ( rawExonBlocks[i].rightSplice != -1 )
						exonBlocks[k].rightSplice = rawExonBlocks[i].rightSplice ;
				}
				else
				{
					exonBlocks.push_back( rawExonBlocks[i] ) ;
					++k ;
				}
			}
		}
} ;

#endif

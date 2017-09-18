// The class manage the blocks
// Li Song

#ifndef _LSONG_CLASSES_BLOCKS_HEADER
#define _LSONG_CLASSES_BLOCKS_HEADER

#include <stdlib.h> 
#include <vector>
#include <map>
#include <assert.h>
#include <math.h>
#include <set>
#include <inttypes.h>

#include "defs.h"

extern bool VERBOSE ;
extern FILE *fpOut ;

extern int gMinDepth ; 


struct _splitSite // means the next position belongs to another block
{
	int64_t pos ;
	char strand ;
	int chrId ;
	int type ; // 1-start of an exon. 2-end of an exon.
	
	int support ;
	int uniqSupport ;

	int64_t oppositePos ; // the position of the other sites to form the intron.
} ;

struct _block
{
	int chrId ;
	int contigId ;
	int64_t start, end ;
	int64_t leftSplice, rightSplice ; // Some information about the left splice site and right splice site.
					  // record the leftmost and rightmost coordinates of a splice site within this block 
				          //or the length of read alignment support the splice sites. Or the number of support.
	int64_t depthSum ;

	int leftType, rightType ; //0-soft boundary, 1-start of an exon, 2-end of an exon.
	
	double leftRatio, rightRatio ; // the adjusted ratio-like value to the left and right anchor subexons.      
	double ratio ;

	char leftStrand, rightStrand ;

	int *depth ;

	int prevCnt, nextCnt ; // Some connection information for the subexons.
	int *prev ;
	int *next ;
} ;

class Blocks
{
	private:
		std::map<int, int> exonBlocksChrIdOffset ;

		int64_t Overlap( int64_t s0, int64_t e0, int64_t s1, int64_t e1, int64_t &s, int64_t &e )
		{
			s = e = -1 ;
			if ( e0 < s1 || s0 > e1 )
				return 0 ;
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
		
		void BuildBlockChrIdOffset()
		{
			// Build the map for the offsets of chr id in the exonBlock list.
			exonBlocksChrIdOffset[ exonBlocks[0].chrId] = 0 ;
			int cnt = exonBlocks.size() ;
			for ( int i = 1 ; i < cnt ; ++i )
			{
				if ( exonBlocks[i].chrId != exonBlocks[i - 1].chrId )
					exonBlocksChrIdOffset[ exonBlocks[i].chrId ] = i ;
			}
		}

		void AdjustAndCreateExonBlocks( int tag, std::vector<struct _block> &newExonBlocks )
		{
			int i, j ;
			if ( exonBlocks[tag].depth != NULL )
			{	
				// Convert the increment and decrement into actual depth.
				int len = exonBlocks[tag].end - exonBlocks[tag].start + 1 ;
				int *depth = exonBlocks[tag].depth ;
				for ( i = 1 ; i < len ; ++i )
					depth[i] = depth[i - 1] + depth[i] ;

				struct _block island ; // the portion created by the hollow.
				island.start = island.end = -1 ;
				island.depthSum = 0 ;

				/*if ( exonBlocks[tag].start == 1562344 )
				{
					for ( i = 0 ; i < len ; ++i )
						printf( "%d\n", depth[i] ) ;
				}*/
				// Adjust boundary accordingly. 
				int64_t adjustStart = exonBlocks[tag].start ; 
				int64_t adjustEnd = exonBlocks[tag].end ;
				
				if ( exonBlocks[tag].leftType == 0 && exonBlocks[tag].rightType != 0 )		
				{
					for ( i = len - 1 ; i >= 0 ; --i )
						if ( depth[i] < gMinDepth )
							break ;
					++i ;
					if ( exonBlocks[tag].rightType == 2 && i + exonBlocks[tag].start > exonBlocks[tag].rightSplice )
						i = exonBlocks[tag].rightSplice - exonBlocks[tag].start ;
					adjustStart = i  + exonBlocks[tag].start ;

					if ( adjustStart > exonBlocks[tag].start )
					{
						int s, e ;
						// firstly [s,e] is the range of depth array.
						for ( s = 0 ; s < i ; ++s )
							if ( depth[s] >= gMinDepth )
								break ;
						for ( e = i - 1 ; e >= s ; -- e )
							if ( depth[e] >= gMinDepth )
								break ;
						if ( e >= s )
						{
							island = exonBlocks[tag] ;
							island.depthSum = 0 ;
							island.leftType = island.rightType = 0 ;
							for ( j = s ; j <= e ; ++j )
								island.depthSum += depth[j] ;
							island.start = s + exonBlocks[tag].start ; // offset the coordinate.
							island.end = e + exonBlocks[tag].start ;
						}
					}
				}
				else if ( exonBlocks[tag].leftType != 0 && exonBlocks[tag].rightType == 0 )
				{
					for ( i = 0 ; i < len ; ++i )	
						if ( depth[i] < gMinDepth )			
							break ;
					--i ;
					if ( exonBlocks[tag].leftType == 1 && i + exonBlocks[tag].start < exonBlocks[tag].leftSplice )
						i = exonBlocks[tag].leftSplice - exonBlocks[tag].start ;
					adjustEnd = i + exonBlocks[tag].start ; 

					if ( adjustEnd < exonBlocks[tag].end )
					{
						int s, e ;
						// firstly [s,e] is the range of depth array.
						for ( s = i + 1 ; s < len ; ++s )
							if ( depth[s] >= gMinDepth )
								break ;
						for ( e = len - 1 ; e >= s ; --e )
							if ( depth[e] >= gMinDepth )
								break ;
						if ( e >= s )
						{
							island = exonBlocks[tag] ;
							island.depthSum = 0 ;
							island.leftType = island.rightType = 0 ;
							for ( j = s ; j <= e ; ++j )
								island.depthSum += depth[j] ;
							island.start = s + exonBlocks[tag].start ; // offset the coordinate.
							island.end = e + exonBlocks[tag].start ;
						}
					}
				}
				else if ( exonBlocks[tag].leftType == 0 && exonBlocks[tag].rightType == 0 )
				{
					for ( i = 0 ; i < len ; ++i )			
						if ( depth[i] >= gMinDepth )
							break ;
					adjustStart = i + exonBlocks[tag].start ;

					for ( i = len - 1 ; i >= 0 ; --i )
						if ( depth[i] >= gMinDepth )
							break ;
					adjustEnd = i + exonBlocks[tag].start ;
				}

				int lostDepthSum = 0 ;
				for ( i = exonBlocks[tag].start ; i < adjustStart ; ++i  )
					lostDepthSum += depth[i - exonBlocks[tag].start ] ;
				for ( i = adjustEnd + 1 ; i < exonBlocks[tag].end ; ++i  )
					lostDepthSum += depth[i - exonBlocks[tag].start ] ;
				exonBlocks[tag].depthSum -= lostDepthSum ;

				delete[] exonBlocks[tag].depth ;
				exonBlocks[tag].depth = NULL ;
				
				//if ( exonBlocks[tag].start == 1562344 )
				//	printf( "%d %d\n", adjustStart, adjustEnd ) ;
				if ( ( len > 1 && adjustEnd - adjustStart + 1 <= 1 ) || ( adjustEnd - adjustStart + 1 <= 0 ) )
					return ;
				exonBlocks[tag].start = adjustStart ;
				exonBlocks[tag].end = adjustEnd ;
				
				if ( island.start != -1 && island.end < exonBlocks[tag].start )
					newExonBlocks.push_back( island ) ;

				newExonBlocks.push_back( exonBlocks[tag] ) ;
				if ( island.start != -1 && island.start > exonBlocks[tag].end )
					newExonBlocks.push_back( island ) ;
			}
		}
	public:
		std::vector<struct _block> exonBlocks ;

		Blocks() { } 	
		~Blocks() 
		{
			int blockCnt = exonBlocks.size() ;
			for ( int i = 0 ; i <  blockCnt ; ++i ) 
			{
				if ( exonBlocks[i].next != NULL )
				{
					delete[] exonBlocks[i].next ;
				}
				if ( exonBlocks[i].prev != NULL )
				{
					delete[] exonBlocks[i].prev ;
				}
			}
		}
		
		double GetAvgDepth( const struct _block &block )
		{
			return block.depthSum / (double)( block.end - block.start + 1 ) ;
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
					//	printf( "hi %d %s %d %d\n", i, alignments.GetReadId(), segments[i].a, segments[i].b ) ;
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
			  printf( "%d %d\n", exonBlocks[i].start, exonBlocks[i].end ) ;
			  }*/

			if ( exonBlocks.size() > 0 )
			{
				BuildBlockChrIdOffset() ;
				int cnt = exonBlocks.size() ;
				for ( int i = 0 ; i < cnt ; ++i )
				{
					exonBlocks[i].contigId = exonBlocks[i].chrId ;
					
					exonBlocks[i].leftType = 0 ;
					exonBlocks[i].rightType = 0 ;
					exonBlocks[i].depth = NULL ;
					exonBlocks[i].nextCnt = 0 ;
					exonBlocks[i].prevCnt = 0 ;
					exonBlocks[i].leftStrand = '.' ;
					exonBlocks[i].rightStrand = '.' ;
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
				//printf( "%s %s: %d %d\n", alignments.GetChromName( splitSites[tag].chrId ), alignments.GetChromName( rawExonBlocks[i].chrId ), tag, l ) ;
				for ( j = tag ; j <= l ; ++j )
				{
					int leftType = 0 ;
					int rightType = 0 ;
					if ( j > tag )	
					{
						start = splitSites[j - 1].pos ; // end is from previous stage
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

					if ( leftType == 2 )
						++start ;
					if ( rightType == 1 )
						--end ;

					struct _block tmpB ;
					tmpB = rawExonBlocks[i] ;
					tmpB.start = start ;
					tmpB.end = end ;
					tmpB.depthSum = 0 ;
					tmpB.ratio = 0 ;
					//printf( "\t%lld %lld\n", start, end ) ;
					tmpB.leftType = leftType ;
					tmpB.rightType = rightType ;
					tmpB.depth = NULL ;
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
					if ( start > end )
						//|| ( tmpB.start == rawExonBlocks[i].start && tmpB.end != rawExonBlocks[i].end && tmpB.end - tmpB.start + 1 <= 7 ) 
						//|| ( tmpB.end == rawExonBlocks[i].end && tmpB.start != rawExonBlocks[i].start && tmpB.end - tmpB.start + 1 <= 7  )) // the case of too short overhang
					{
						exonBlocks.pop_back() ;
					}
					/*else if ( start == end )
					{

					}*/
				}
			}
			BuildBlockChrIdOffset() ;
		}

		void ComputeDepth( Alignments &alignments ) 
		{
			// Go through the alignment list again to fill in the depthSum;
			int i ;
			int j ;
			int tag = 0 ;
			int blockCnt = exonBlocks.size() ;

			std::vector<struct _block> newExonBlocks ;
			while ( alignments.Next() )
			{
				int segCnt = alignments.segCnt ;
				struct _pair *segments = alignments.segments ;

				while ( tag < blockCnt && ( exonBlocks[tag].chrId < alignments.GetChromId() || 
							( exonBlocks[tag].chrId == alignments.GetChromId() && exonBlocks[tag].end < segments[0].a ) ) )
				{
					AdjustAndCreateExonBlocks( tag, newExonBlocks ) ;
					++tag ;
				}
				for ( i = 0 ; i < segCnt ; ++i )
				{
					for ( j = tag ; j < blockCnt && ( exonBlocks[j].chrId == alignments.GetChromId() && exonBlocks[j].start <= segments[i].b ) ; ++j )
					{
						int64_t s = -1, e = -1 ;
						exonBlocks[j].depthSum += Overlap( segments[i].a, segments[i].b, exonBlocks[j].start, exonBlocks[j].end, s, e ) ; 
						if ( s == -1 )
							continue ;

						if ( exonBlocks[j].depth == NULL )
						{
							int len = exonBlocks[j].end - exonBlocks[j].start + 1 ;
							exonBlocks[j].depth = new int[len + 1] ;
							memset( exonBlocks[j].depth, 0, sizeof( int ) * ( len + 1 ) ) ;
							exonBlocks[j].leftSplice = exonBlocks[j].start ;
							exonBlocks[j].rightSplice = exonBlocks[j].end ;
						}
						int *depth = exonBlocks[j].depth ;
						++depth[s - exonBlocks[j].start] ;
						--depth[e - exonBlocks[j].start + 1] ; // notice the +1 here, since the decrease of coverage actually happens at next base.

						// record the longest alignment stretch support the splice site.
						if ( exonBlocks[j].leftType == 1 && segments[i].a == exonBlocks[j].start 
							&& e > exonBlocks[j].leftSplice )
						{
							exonBlocks[j].leftSplice = e ;	
						}
						if ( exonBlocks[j].rightType == 2 && segments[i].b == exonBlocks[j].end 
							&& s < exonBlocks[j].rightSplice )
						{
							exonBlocks[j].rightSplice = s ;	
						}
					}
				}
			}

			for ( ; tag < blockCnt ; ++tag )
				AdjustAndCreateExonBlocks( tag, newExonBlocks ) ;
			exonBlocks.clear() ;

			// Due to multi-alignment, we may skip some alignments that determines leftSplice and
			// rightSplice, hence filtered some subexons. As a result, some subexons' anchor subexon will disapper
			// and we need to change its boundary type.
			blockCnt = newExonBlocks.size() ;
			for ( i = 0 ; i < blockCnt ; ++i )
			{
				if ( ( i == 0 && newExonBlocks[i].leftType == 2 ) ||
					( i > 0 && newExonBlocks[i].leftType == 2 && newExonBlocks[i - 1].end + 1 != newExonBlocks[i].start ) )
				{
					newExonBlocks[i].leftType = 0 ;
				}

				if ( ( i == blockCnt - 1 && newExonBlocks[i].rightType == 1 ) ||
					( i < blockCnt - 1 && newExonBlocks[i].rightType == 1 && 
						newExonBlocks[i].end + 1 != newExonBlocks[i + 1].start ) )
					newExonBlocks[i].rightType = 0 ;
			}
			exonBlocks = newExonBlocks ;
		}

		// If two blocks whose soft boundary are close to each other, we can merge them.
		void MergeNearBlocks()
		{
			std::vector<struct _block> rawExonBlocks = exonBlocks ;
			int i, k ;
			int bsize = rawExonBlocks.size() ; 
			int *newIdx = new int[bsize] ; // used for adjust prev and next. Note that by merging, it won't change the number of prevCnt or nextCnt.
			
			exonBlocks.clear() ;
			exonBlocks.push_back( rawExonBlocks[0] ) ;
			k = 0 ;
			newIdx[0] = 0 ;
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
				newIdx[i] = k ;
			}
			BuildBlockChrIdOffset() ;
			bsize = exonBlocks.size() ;
			for ( i = 0 ; i < bsize ; ++i )
			{
				int cnt = exonBlocks[i].prevCnt ;
				for ( k = 0 ; k < cnt ; ++k )
					exonBlocks[i].prev[k] = newIdx[ exonBlocks[i].prev[k] ] ;

				cnt = exonBlocks[i].nextCnt ;
				for ( k = 0 ; k < cnt ; ++k )
					exonBlocks[i].next[k] = newIdx[ exonBlocks[i].next[k] ] ;
			}
			delete[] newIdx ;
		}

		void AddIntronInformation( std::vector<struct _splitSite> &sites )
		{
			// Add the connection information, support information and the strand information.
			int i, j, k, tag ;
			tag = 0 ;
			int scnt = sites.size() ;
			int exonBlockCnt = exonBlocks.size() ;
			bool flag ;

			for ( i = 0 ; i < exonBlockCnt ; ++i )
			{
				exonBlocks[i].prevCnt = exonBlocks[i].nextCnt = 0 ;
				exonBlocks[i].prev = exonBlocks[i].next = NULL ;

				exonBlocks[i].leftStrand = exonBlocks[i].rightStrand = '.' ;

				exonBlocks[i].leftSplice = exonBlocks[i].rightSplice = 0 ;
			}
			
			// Now, we add the region id
			int regionId = 0 ;
			for ( i = 0 ; i < exonBlockCnt ; )
			{
				for ( j = i + 1 ; j < exonBlockCnt ; ++j )
					if ( exonBlocks[j].start > exonBlocks[j - 1].end + 1 )
						break ;
				for ( k = i ; k < j ; ++k )
					exonBlocks[k].contigId = regionId ;
				++regionId ;
				i = j ;
			}

			for ( i = 0 ; i < scnt ; )	
			{
				for ( j = i ; j < scnt ; ++j )		
				{
					if ( sites[j].chrId != sites[i].chrId ||
							sites[j].pos != sites[i].pos ||
							sites[j].type != sites[i].type )
						break ;
				}
				// [i,j-1] are the indices of the sites with same coordinate

				int cnt = j - i ;
				// Locate the first subexon that can overlap with this site
				while ( tag < exonBlockCnt )
				{
					if ( ( exonBlocks[tag].chrId == sites[i].chrId && exonBlocks[tag].end >= sites[i].pos ) 
						|| exonBlocks[tag].chrId > sites[i].chrId )
						break ;
					++tag ;
				}
				flag = false ;
				for ( k = tag; k < exonBlockCnt ; ++k )
				{
					if ( exonBlocks[tag].start > sites[i].pos ||
						exonBlocks[tag].chrId > sites[i].chrId )
						break ;

					if ( exonBlocks[tag].start == sites[i].pos ||
						exonBlocks[tag].end == sites[i].pos )
					{
						flag = true ;
						break ;
					}
				}
				
				if ( !flag )
				{
					i = j ; 
					continue ;
				}

				if ( sites[i].type == 1 && exonBlocks[tag].start == sites[i].pos )
				{
					exonBlocks[tag].prevCnt = 0 ;
					exonBlocks[tag].prev = new int[cnt] ;
					exonBlocks[tag].leftStrand = sites[i].strand ;

					// And we also need to put the "next" here.
					// Here we assume the oppositePos sorted in increasing order
					k = tag - 1 ;
					for ( int l = j - 1 ; l >= i ; --l )
					{
						for ( ; k >= 0 ; --k )			
						{
							if ( exonBlocks[k].end < sites[l].oppositePos || exonBlocks[k].chrId != sites[l].chrId )
								break ;
							if ( exonBlocks[k].end == sites[l].oppositePos )
							{
								exonBlocks[k].next[ exonBlocks[k].nextCnt ] = tag ;
								exonBlocks[tag].prev[ exonBlocks[tag].prevCnt] = k ; // Note that the prev are sorted in decreasing order
								++exonBlocks[k].nextCnt ;
								++exonBlocks[tag].prevCnt ;
								
								if ( exonBlocks[k].contigId != exonBlocks[tag].contigId )
								{
									exonBlocks[tag].leftSplice += sites[l].support ;
									exonBlocks[k].rightSplice += sites[l].support ;
								}
								break ;
							}
						}
					}
				}
				else if ( sites[i].type == 2 && exonBlocks[tag].end == sites[i].pos )
				{
					exonBlocks[tag].nextCnt = 0 ; // cnt ; it should reach cnt after putting the ids in
					exonBlocks[tag].next = new int[cnt] ;
					exonBlocks[tag].rightStrand = sites[i].strand ;
				}

				i = j ;
			}
			// Reverse the prev list
			for ( i = 0 ; i < exonBlockCnt ; ++i )
			{
				for ( j = 0, k = exonBlocks[i].prevCnt - 1 ; j < k ; ++j, --k )
				{
					int tmp = exonBlocks[i].prev[j] ;
					exonBlocks[i].prev[j] = exonBlocks[i].prev[k] ;
					exonBlocks[i].prev[k] = tmp ;
				}
			}

			// If all of the subexon anchored the intron are filtered, we need to change the type of the other anchor.
			for ( i = 0 ; i < exonBlockCnt ; ++i )		
			{
				if ( exonBlocks[i].leftType == 1 && exonBlocks[i].prevCnt == 0 )	
				{
					exonBlocks[i].leftType = 0 ;
					if ( i > 0 && exonBlocks[i - 1].rightType == 1 )	
						exonBlocks[i - 1].rightType = 0 ;
				}
				if ( exonBlocks[i].rightType == 2 && exonBlocks[i].nextCnt == 0 )
				{
					exonBlocks[i].rightType = 0 ;
					if ( i < exonBlockCnt - 1 && exonBlocks[i + 1].leftType == 2 )
						exonBlocks[i + 1].leftType = 0 ;
				}
			}
			MergeNearBlocks() ;
		}

		double PickLeftAndRightRatio( double l, double r )
		{
			if ( l >= 0 && r >= 0 )
				return l < r ? l : r ;
			else if ( l < 0 && r < 0 )
				return -1 ;
			else if ( l < 0 )
				return r ;
			else
				return l ;
		}
			
		double PickLeftAndRightRatio( const struct _block &b )
		{
			return PickLeftAndRightRatio( b.leftRatio, b.rightRatio ) ;
		}

		void ComputeRatios()
		{
			int i, j ;
			int exonBlockCnt = exonBlocks.size() ;		
			for ( i = 0 ; i < exonBlockCnt ; ++i )
			{
				exonBlocks[i].leftRatio = -1 ;
				exonBlocks[i].rightRatio = -1 ;
				if ( exonBlocks[i].leftType == 2 && exonBlocks[i].rightType == 1 )
				{
					// ]...[
					double anchorAvg = 0 ;
					double avgDepth = GetAvgDepth( exonBlocks[i] ) ;
					if ( i >= 1 && exonBlocks[i - 1].chrId == exonBlocks[i].chrId )
					{
						anchorAvg = GetAvgDepth( exonBlocks[i - 1] ) ;	
						if ( anchorAvg > 1 && avgDepth > 1 )
							exonBlocks[i].leftRatio = ( avgDepth - 1 ) / ( anchorAvg - 1 ) ;  
					}
					if ( i < exonBlockCnt - 1 && exonBlocks[i + 1].chrId == exonBlocks[i].chrId )
					{
						anchorAvg = GetAvgDepth( exonBlocks[i + 1] ) ;	
						if ( anchorAvg > 1 && avgDepth > 1 )
							exonBlocks[i].rightRatio = ( avgDepth - 1 ) / ( anchorAvg - 1 ) ;  
					}
				}

				if ( exonBlocks[i].leftType == 0 && exonBlocks[i].rightType == 1 ) 
				{
					double avgDepth = GetAvgDepth( exonBlocks[i] ) ;
					double anchorAvg = GetAvgDepth( exonBlocks[i + 1] ) ;
					if ( avgDepth > 1 && anchorAvg > 1 )
						exonBlocks[i].rightRatio = ( avgDepth - 1 ) / ( anchorAvg - 1 ) ;
				}

				/*if ( exonBlocks[i].leftType == 1 )
				{
					// For the case (...[, the leftRatio is actuall the leftratio of the subexon on its right. 	
					int len = 0 ;
					double depthSum = 0 ;
					int tag = i ;
					
					for ( j = 0 ; j < exonBlocks[tag].prevCnt ; ++j )
					{
						int k = exonBlocks[tag].prev[j] ;
						len += ( exonBlocks[k].end - exonBlocks[k].start + 1 ) ;		
						depthSum += exonBlocks[k].depthSum ;	
					}
					double otherAvgDepth = depthSum / len ;
					double avgDepth = GetAvgDepth( exonBlocks[tag] ) ;

					// adjust avgDepth by looking into the following exon blocks(subexon) in the same region whose left side is not hard end
					for ( j = tag + 1 ; j < exonBlockCnt ; ++j )
					{
						if ( exonBlocks[j].start != exonBlocks[j - 1].end + 1 || exonBlocks[j].leftType == 1 )
							break ;
						double tmp = GetAvgDepth( exonBlocks[j] ) ;
						if ( tmp > avgDepth )
							avgDepth = tmp ;
					}

					if ( avgDepth < 1 )
						avgDepth = 1 ;
					if ( otherAvgDepth < 1 )
						otherAvgDepth = 1 ;

					//exonBlocks[i].leftRatio = pow( log( avgDepth ) / log(2.0), 2.0 / 3.0 ) - pow( log( otherAvgDepth ) / log( 2.0 ), 2.0 / 3.0 );
					exonBlocks[i].leftRatio = sqrt( avgDepth ) - log( avgDepth ) - ( sqrt( otherAvgDepth ) + log( otherAvgDepth ) ) ;

				}*/
				if ( exonBlocks[i].rightType == 0 && exonBlocks[i].leftType == 2  ) 
				{
					// for overhang subexon.
					double avgDepth = GetAvgDepth( exonBlocks[i] ) ;
					double anchorAvg = GetAvgDepth( exonBlocks[i - 1] ) ;
					if ( avgDepth > 1 && anchorAvg > 1 )
						exonBlocks[i].leftRatio = ( avgDepth - 1 ) / ( anchorAvg - 1 ) ;

				}
				/*if ( exonBlocks[i].rightType == 2 ) 
				{
					int len = 0 ;
					double depthSum = 0 ;
					int tag = i ;
					
					for ( j = 0 ; j < exonBlocks[tag].nextCnt ; ++j )
					{
						int k = exonBlocks[tag].next[j] ;
						len += ( exonBlocks[k].end - exonBlocks[k].start + 1 ) ;		
						depthSum += exonBlocks[k].depthSum ;	
					}
					double otherAvgDepth = depthSum / len ;
					double avgDepth = GetAvgDepth( exonBlocks[tag] ) ;
					
					// adjust avgDepth by looking into the earlier exon blocks(subexon) in the same region whose right side is not hard end
					for ( j = tag - 1 ; j >= 0 ; --j )
					{
						if ( exonBlocks[j].end + 1 != exonBlocks[j + 1].start || exonBlocks[j].rightType == 2 )
							break ;
						double tmp = GetAvgDepth( exonBlocks[j] ) ;
						if ( tmp > avgDepth )
							avgDepth = tmp ;
					}

					if ( avgDepth < 1 )
						avgDepth = 1 ;
					if ( otherAvgDepth < 1 )
						otherAvgDepth = 1 ;

					exonBlocks[i].rightRatio = sqrt( avgDepth ) - log( avgDepth ) - ( sqrt( otherAvgDepth ) + log( otherAvgDepth ) ) ;
					//pow( log( avgDepth ) / log(2.0), 2.0 / 3.0 ) - pow( log( otherAvgDepth ) / log(2.0), 2.0 / 3.0 );
				}*/
				// The remaining case the islands, (...)
			}

			// go through each region of subexons, and only compute the ratio for the first and last hard boundary
			for ( i = 0 ; i < exonBlockCnt ; )
			{
				// the blocks from [i,j) corresponds to a region, namely consecutive subexons.
				for ( j = i + 1 ; j < exonBlockCnt ; ++j )
					if ( exonBlocks[j].contigId != exonBlocks[i].contigId )
						break ;

				int leftSupport = 0 ;
				int rightSupport = 0 ;
				int leftTag = j, rightTag = i - 1  ;
				for ( int k = i ; k < j ; ++k )
				{
					if ( exonBlocks[k].leftType == 1 && exonBlocks[k].leftSplice > 0 )	
					{
						leftSupport += exonBlocks[k].leftSplice ;
						if ( k < leftTag )
							leftTag = k ;
					}
					if ( exonBlocks[k].rightType == 2 && exonBlocks[k].rightSplice > 0 )
					{
						rightSupport += exonBlocks[k].rightSplice ;
						if ( k > rightTag )
							rightTag = k ;
					}
				}

				if ( leftSupport > 0 && rightSupport > 0 )
				{
					// when the right splice support is much greater than that of the left side, there should be a soft boundary for the left side.
					// Wether we want to include this soft boundary or not will be determined in class, when it looked at whether the overhang block
					// 	should be kept or not.
					exonBlocks[ leftTag ].leftRatio = sqrt( rightSupport ) - log( rightSupport ) - ( sqrt( leftSupport ) + log( leftSupport ) ) ;
					exonBlocks[ rightTag ].rightRatio = sqrt( leftSupport ) - log( leftSupport ) - ( sqrt( rightSupport ) + log( rightSupport ) ) ;
				}

				i = j ; // don't forget this.
			}
		}
} ;

#endif

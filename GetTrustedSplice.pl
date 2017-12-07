#!/bin/perl

use strict ;
use warnings ;

die "usage: a.pl path_to_list_of_splice_file > trusted.splice\n" if ( @ARGV == 0 ) ;

my %spliceSupport ;
my %spliceUniqSupport ;
my %spliceSecSupport ;
my %uniqSpliceSites ;

open FP1, $ARGV[0] ;
my $sampleCnt = 0 ;
while ( <FP1> )
{
	chomp ;
	open FP2, $_ ;
	++$sampleCnt ;
	while ( <FP2> )
	{
		chomp ;
		my $line = $_ ;
		my @cols = split /\s+/, $line ;
		my $key = $cols[0]." ".$cols[1]." ".$cols[2]." ".$cols[4] ;
		if ( $cols[3] <= 0 )
		{
			$cols[3] = 0.1 ;
		}

		if ( ! defined $spliceSupport{$key} )
		{
			$spliceSupport{ $key } = $cols[3] ;
			$spliceUniqSupport{ $key } = $cols[5] ;
			$spliceSecSupport{ $key } = $cols[6] ;
		}
		else
		{
			$spliceSupport{ $key } += $cols[3] ;
			$spliceUniqSupport{ $key } += $cols[5] ;
			$spliceSecSupport{ $key } += $cols[6] ;
		}
		
		for ( my $i = 1 ; $i <=2 ; ++$i )
		{
			$key = $cols[0]." ".$cols[$i] ;
			if ( defined $uniqSpliceSites{ $key } && $uniqSpliceSites{ $key } != $cols[2 - $i + 1] )
			{
				$uniqSpliceSites{ $key } = -1 ;
			}
			else
			{
				$uniqSpliceSites{ $key } = $cols[2 - $i + 1] ;
			}
		}
	}
	close FP2 ;
}
close FP1 ;

foreach my $key (keys %spliceSupport)
{
	next if ( $spliceSupport{ $key } / $sampleCnt < 0.5 ) ;
	next if ( $spliceUniqSupport{$key} / ( $spliceSecSupport{$key} + $spliceUniqSupport{$key} ) < 0.01 ) ;
	my @cols = split /\s+/, $key ;
	my $flag = 0 ;
	if ( $cols[2] - $cols[1] + 1 > 10000 )
	{
		$flag = 1 if ( $spliceSupport{ $key } / $sampleCnt < 1 ) ;
	}
	if ( $cols[2] - $cols[1] + 1 > 100000 )
	{
		$flag = 1 if ( $spliceUniqSupport{$key} / ( $spliceSecSupport{$key} + $spliceUniqSupport{$key} ) < 0.1 ) ;
	}
	if ( $flag == 1 && ( $uniqSpliceSites{ $cols[0]." ".$cols[1] } == -1 || $uniqSpliceSites{ $cols[0]." ".$cols[2] } == -1 ) )
	{
		next ;
	}
	print $cols[0], " ", $cols[1], " ", $cols[2], " 10 ", $cols[3], " 10 0 0 0\n" ;
}

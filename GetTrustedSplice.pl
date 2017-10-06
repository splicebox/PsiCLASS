#!/bin/perl

use strict ;
use warnings ;

die "usage: a.pl path_to_list_of_splice_file > trusted.splice\n" if ( @ARGV == 0 ) ;

my %spliceSupport ;

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
			$cols[3] = 1 ;
		}

		if ( ! defined $spliceSupport{$key} )
		{
			$spliceSupport{ $key } = $cols[3] ;
		}
		else
		{
			$spliceSupport{ $key } += $cols[3] ;
		}
	}
	close FP2 ;
}
close FP1 ;

foreach my $key (keys %spliceSupport)
{
	next if ( $spliceSupport{ $key } / $sampleCnt < 0.5 ) ;
	my @cols = split /\s+/, $key ;
	print $cols[0], " ", $cols[1], " ", $cols[2], " 10 ", $cols[3], " 10 0 0 0\n" ;
}

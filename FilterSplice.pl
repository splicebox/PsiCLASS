#!/bin/perl

use strict ;
use warnings ;

die "usage: a.pl a.raw_splice trusted.plisce > filtered.splice\n" if ( @ARGV == 0 ) ;

my %trustedSplices ;

open FP1, $ARGV[1] ;
while ( <FP1> )
{
	chomp ;
	my @cols = split /\s+/, $_ ;
	my $key = $cols[0]." ".$cols[1]." ".$cols[2] ;
	$trustedSplices{ $key } = 1 ;
}
close FP1 ;

open FP1, $ARGV[0] ;
while ( <FP1> )
{
	chomp ;
	my @cols = split /\s+/, $_ ;
	my $key = $cols[0]." ".$cols[1]." ".$cols[2] ;
	next if ( !defined $trustedSplices{ $key } ) ;
	
	if ( $cols[3] <= 0 )
	{
		print $cols[0], " ", $cols[1], " ", $cols[2], " 1 ", $cols[4], " 1 0 0 0\n" ;
	}
	else
	{
		print $_, "\n" ;
	}

}
close FP1 ;


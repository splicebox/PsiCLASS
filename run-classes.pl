#!/usr/bin/env perl

use strict ;
use warnings ;

use Cwd 'cwd' ;
use Cwd 'abs_path' ;
use File::Basename;

die "Usage: perl run-classes.pl [OPTIONS]\n".
    "Required:\n".
    "\t--lb STRING: the path to the file containing the alignments bam files\n".
    "Optional:\n".
    "\t-s STRING: the path to the trusted splice sites file (default: not used)\n".
    "\t-p STRING: the prefix of output files (default: classes)\n". 
    "\t--stage NUM:  (default: 0)\n".
    "\t\t0-start from beginning (building splice sites for each sample)\n".
    "\t\t1-start from building subexon files for each sample\n".
    "\t\t2-start from combining subexon files\n".
    "\t\t3-start from assembly the transcripts\n"
    if ( @ARGV == 0 ) ;

my $WD = dirname( abs_path( $0 ) ) ;

my $i ;
my $cmd ;
my $prefix = "classes_" ;

sub system_call
{
	print $_[0], "\n" ;
	system( $_[0] ) ;
}

# Process the arguments
my @bamFiles ;
my $spliceFile = "" ;
my $bamFileList ;
my $stage = 0 ;
for ( $i = 0 ; $i < @ARGV ; ++$i )
{
	if ( $ARGV[$i] eq "--lb" )
	{
		$bamFileList = $ARGV[$i + 1] ;
		open FP1, $ARGV[$i + 1] ;	
		while ( <FP1> )
		{
			chomp ;
			push @bamFiles, $_ ;
		}
		++$i ;
	}
	elsif ( $ARGV[$i] eq "-s" )
	{
		$spliceFile = $ARGV[$i + 1] ;		
		++$i ;
	}
	elsif ( $ARGV[$i] eq "-p" )
	{
		$prefix = $ARGV[$i + 1] ;
		if ( substr( $prefix, -1 ) ne "_" )
		{
			$prefix .= "_" ;
		}
		++$i ;
	}
	elsif ( $ARGV[$i] eq "--stage" )
	{
		$stage = $ARGV[$i + 1] ;
		++$i ;
	}
}
if ( scalar( @bamFiles ) == 0 )
{
	die "Must use option --lb to specify the list of bam files.\n" ;
}

# Generate the splice file for each bam file.
if ( $stage <= 0 )
{
	system_call( "echo -n > ${prefix}splice.list" ) ;
	for ( $i = 0 ; $i < @bamFiles ; ++$i )
	{
		system_call( "$WD/junc ".$bamFiles[$i]." -a > ${prefix}bam_$i.raw_splice" ) ;
		#if ( $spliceFile ne "" )
		#{
		#	system_call( "perl $WD/ManipulateIntronFile.pl $spliceFile ${prefix}bam_$i.raw_splice > ${prefix}bam_$i.splice" ) ;
		#}
		#else
		#{
#system_call( "awk \'{if (\$6>1) print;}\' ${prefix}bam_$i.raw_splice > ${prefix}bam_$i.splice" ) ;
		#	system_call( "mv ${prefix}bam_$i.raw_splice ${prefix}bam_$i.splice" ) ;
		#}

		system_call( "echo ${prefix}bam_$i.raw_splice >> ${prefix}splice.list" )
	}
	
	if ( $spliceFile ne "" )
	{
		for ( $i = 0 ; $i < @bamFiles ; ++$i )
		{
			system_call( "perl $WD/FilterSplice.pl ${prefix}bam_$i.raw_splice $spliceFile > ${prefix}bam_$i.splice" ) ;
		}
	}
	else
	{
		system_call( "perl $WD/GetTrustedSplice.pl ${prefix}splice.list > ${prefix}bam.trusted_splice" ) ;
		for ( $i = 0 ; $i < @bamFiles ; ++$i )
		{
			system_call( "perl $WD/FilterSplice.pl ${prefix}bam_$i.raw_splice ${prefix}bam.trusted_splice > ${prefix}bam_$i.splice" ) ;
		}
	}
}

# Get subexons from each bam file
if ( $stage <= 1 )
{
	open FPls, ">${prefix}subexon.list" ;
	for ( $i = 0 ; $i < @bamFiles ; ++$i )
	{
		system_call( "$WD/subexon-info ".$bamFiles[$i]." ${prefix}bam_$i.splice > ${prefix}subexon_$i.out" ) ;	
		print FPls "${prefix}subexon_$i.out\n" ;
	}
}

# combine the subexons.
if ( $stage <= 2 )
{
	$cmd = "$WD/combine-subexons --ls ${prefix}subexon.list > ${prefix}subexon_combined.out" ;
	system_call( "$cmd" ) ;
}

# Run classes
if ( $stage <= 3 )
{
	my $trimPrefix = substr( $prefix, 0, -1 ) ;
	$cmd = "$WD/classes --lb $bamFileList -s ${prefix}subexon_combined.out -o ${trimPrefix} > classes.log" ;
	system_call( "$cmd" ) ;
}

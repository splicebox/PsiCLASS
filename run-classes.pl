#!/usr/bin/env perl

use strict ;
use warnings ;

use Cwd 'cwd' ;
use Cwd 'abs_path' ;
use File::Basename;
use threads ;
use threads::shared ;

die "Usage: perl run-classes.pl [OPTIONS]\n".
    "Required:\n".
    "\t-b STRING: the path to the BAM files. Use comma to separate ultiple bam files\n".
    "\t\tor\n".
    "\t--lb STRING: the path to the file containing the alignments bam files\n".
    "Optional:\n".
    "\t-s STRING: the path to the trusted splice sites file (default: not used)\n".
    "\t-o STRING: the prefix of output files (default: classes)\n". 
    "\t-t INT: number of threads (default: 1)\n".
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
my $numThreads = 1 ;

sub system_call
{
	print $_[0], "\n" ;
	system( $_[0] ) ;
}

# Process the arguments
my @bamFiles : shared ;
my $spliceFile = "" ;
my $bamFileList = "" ;
my $stage = 0 ;
my $classesOpt = "" ;
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
	elsif ( $ARGV[$i] eq "-b" )
	{
		@bamFiles = split /,/, $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "-s" )
	{
		$spliceFile = $ARGV[$i + 1] ;		
		++$i ;
	}
	elsif ( $ARGV[$i] eq "-o" )
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
	elsif ( $ARGV[ $i ] eq "-t" )
	{
		$numThreads = $ARGV[$i + 1] ;
		$classesOpt .= " -t $numThreads" ;
		++$i ;
	}
	else
	{
		die "Unknown argument: ", $ARGV[$i], "\n" ;
	}
}
if ( scalar( @bamFiles ) == 0 )
{
	die "Must use option --lb to specify the list of bam files.\n" ;
}

my $threadLock : shared ;
my @sharedFiles : shared ;
my @threads ;
for ( $i = 0 ; $i < $numThreads ; ++$i )
{
	push @threads, $i ;
}


sub threadRunSplice
{
	my $tid = threads->tid() - 1 ;
	my $i ;
	for ( $i = 0 ; $i < scalar( @bamFiles ) ; ++$i )	
	{
		next if ( ( $i % $numThreads ) != $tid ) ;
		system_call( "$WD/junc ".$bamFiles[$i]." -a > ${prefix}bam_$i.raw_splice" ) ;
	}
}

# Generate the splice file for each bam file.
if ( $stage <= 0 )
{
	if ( $numThreads == 1 )
	{
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
		}
	}
	else
	{
		foreach ( @threads )
		{
			$_ = threads->create( \&threadRunSplice ) ;
		}
		foreach ( @threads )
		{
			$_->join() ;
		}
	}
	
	open FPls, ">${prefix}splice.list" ;
	for ( $i = 0 ; $i < @bamFiles ; ++$i )
	{
		print FPls  "${prefix}bam_$i.raw_splice\n" ;
	}
	close FPls ;

	if ( $spliceFile ne "" )
	{
		for ( $i = 0 ; $i < @bamFiles ; ++$i )
		{
			system_call( "perl $WD/FilterSplice.pl ${prefix}bam_$i.raw_splice $spliceFile > ${prefix}bam_$i.splice" ) ;
		}
	}
	else
	{
		system_call( "$WD/trust-splice ${prefix}splice.list ". $bamFiles[0] ." > ${prefix}bam.trusted_splice" ) ;
		for ( $i = 0 ; $i < @bamFiles ; ++$i )
		{
			system_call( "perl $WD/FilterSplice.pl ${prefix}bam_$i.raw_splice ${prefix}bam.trusted_splice > ${prefix}bam_$i.splice" ) ;
		}
	}
}


# Get subexons from each bam file
sub threadRunSubexonInfo
{
	my $tid = threads->tid() - $numThreads - 1 ;
	my $i ;
	for ( $i = 0 ; $i < scalar( @bamFiles ) ; ++$i )	
	{
		next if ( ( $i % $numThreads ) != $tid ) ;
		system_call( "$WD/subexon-info ".$bamFiles[$i]." ${prefix}bam_$i.splice > ${prefix}subexon_$i.out" ) ;	
	}
}


if ( $stage <= 1 )
{

	if ( $numThreads == 1 )
	{
		for ( $i = 0 ; $i < @bamFiles ; ++$i )
		{
			system_call( "$WD/subexon-info ".$bamFiles[$i]." ${prefix}bam_$i.splice > ${prefix}subexon_$i.out" ) ;	
		}
	}
	else
	{
		print "hi ", scalar( @threads ), " ", scalar( @bamFiles), "\n" ;
		foreach ( @threads )
		{
			$_ = threads->create( \&threadRunSubexonInfo ) ;
		}
		foreach ( @threads )
		{
			$_->join() ;
		}
	}
	
	open FPls, ">${prefix}subexon.list" ;
	for ( $i = 0 ; $i < @bamFiles ; ++$i )
	{
		print FPls "${prefix}subexon_$i.out\n" ;
	}
	close FPls ;
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
	my $bamPath = "" ;
	if ( $bamFileList ne "" )
	{
		$bamPath = " --lb $bamFileList " ;
	}
	else
	{
		foreach my $b (@bamFiles)
		{
			$bamPath .= " -b $b " ;
		}
	}
	$cmd = "$WD/classes $bamPath -s ${prefix}subexon_combined.out -o ${trimPrefix} > classes.log" ;
	system_call( "$cmd" ) ;
}

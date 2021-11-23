#!/usr/bin/perl

use warnings;
use strict;
use diagnostics;
use IO::File;

my $groupingsFile = $ARGV[0];
my $fastaInput    = $ARGV[1];
my $backgroundBF  = $ARGV[2];
my $fh            = new IO::File( $groupingsFile, "r" );

my $line = $fh->getline();

my %sampleToGroup;

#extract sampleID->groupID
while ($line) {
	chomp($line);
	my @tempArr          = split( /\t/, $line );
	my $sampleName       = $tempArr[0];
	my $tempFilename     = "temp.fa";
	my $tempOutputPrefix = "temp";

#GRCh38#0#chr6:31984662-32023790	117.912063381843	-11.4816618982337	14035.0832509091	1
	my $groupID = $tempArr[4];

	#extract sampleID->groupID
	$sampleToGroup{$sampleName} = $groupID;

	#extract each fasta file seperately from input fasta
	my $cmd =
	  'grep -A1 ' . $sampleName . ' ' . $fastaInput . ' > ' . $tempFilename;

	system($cmd);

	#run dwgsim
	my $dwgsimCMD =
	  'dwgsim -1 150 -2 150 -C 100 ' . $tempFilename . ' ' . $tempOutputPrefix;

	system($dwgsimCMD);

#run classification code on each
#/home/cjustin/git/KmerClassifier/src/seqgc -p test -g ../../../C4.grouping.tsv -b ../C4_sgc_bg.bf -i C4_HG00735.2_reads.bfast.fastq.gz -c 0.001 ../C4_fixed.fa
	my $seqgcCMD =
		'seqgc -p '
	  . $tempOutputPrefix . ' -g '
	  . $groupingsFile . ' -b '
	  . $backgroundBF . ' -i '
	  . $tempOutputPrefix
	  . '.bfast.fastq.gz '
	  . $fastaInput;

	my $results        = `$seqgcCMD`;
	my @resultsPerLine = split( /\n/, $results );


	#parse results
	foreach my $resultsLine (@resultsPerLine) {
		#sampleID groupID KLValue correctGroupID
		my @resultsTempArr = split( /\t/, $resultsLine );
		print $sampleName . "\t"
		  . $resultsTempArr[0] . "\t"
		  . $resultsTempArr[1] . "\t"
		  . $groupID . "\n";
	}

	#sampleID bestKL secondBestKL TODO
	
	
	
	$line = $fh->getline();
}
$fh->close();

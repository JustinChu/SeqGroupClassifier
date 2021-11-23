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
my @sampleIDs;

#extract sampleID->groupID
while ($line) {
	chomp($line);
	my @tempArr    = split( /\t/, $line );
	my $sampleName = $tempArr[0];

#GRCh38#0#chr6:31984662-32023790	117.912063381843	-11.4816618982337	14035.0832509091	1
	my $groupID = $tempArr[4];

	#extract sampleID->groupID
	$sampleToGroup{$sampleName} = $groupID;
	push( @sampleIDs, $sampleName );
	$line = $fh->getline();
}
$fh->close();
for ( my $i = 0 ; $i < scalar(@sampleIDs) ; ++$i ) {
	for ( my $j = $i ; $j < scalar(@sampleIDs) ; ++$j ) {
		my $tempFilename     = $sampleIDs[$i] . "_" . $sampleIDs[$j].".fa";
		my $tempOutputPrefix = $sampleIDs[$i] . "_" . $sampleIDs[$j];
		my $sampleID1        = $sampleIDs[$i];
		my $sampleID2        = $sampleIDs[$j];

		#for each combination of samples

		#extract each fasta file seperately from input fasta
		my $cmd1 =
		  'grep -A1 ' . $sampleID1 . ' ' . $fastaInput . ' > ' . $tempFilename;

		system($cmd1);

		my $cmd2 =
		  'grep -A1 ' . $sampleID2 . ' ' . $fastaInput . ' >> ' . $tempFilename;

		system($cmd2);

#run classification code on each
#/home/cjustin/git/KmerClassifier/src/seqgc -p test -g ../../../C4.grouping.tsv -b ../C4_sgc_bg.bf -i C4_HG00735.2_reads.bfast.fastq.gz -c 0.001 ../C4_fixed.fa
		my $seqgcCMD =
			'seqgc --debug -p '
		  . $tempOutputPrefix . ' -g '
		  . $groupingsFile . ' -b '
		  . $backgroundBF . ' -i '
		  . $tempFilename . ' '
		  . $fastaInput;

		
		print STDERR $sampleIDs[$i] . "\t". $sampleIDs[$j]  . "\n" . $seqgcCMD . "\n";
		
		my $results = `$seqgcCMD`;

		my @resultsPerLine = split( /\n/, $results );

		#parse results
		foreach my $resultsLine (@resultsPerLine) {

#0.69149	HG02886#1#h1tg000006l:23458349-23491138	HG03516#2#h2tg000003l:32100047-32132824	7	7
#sampleID1 sampleID2 trueSampleID1 trueSampleID2 groupID1 groupID2 trueGroupID1 trueGroupID2 KLValue
			my @resultsTempArr = split( /\t/, $resultsLine );
			print $resultsTempArr[0] . "\t"
			  . $resultsTempArr[1] . "\t"
			  . $sampleID1 . "\t"
			  . $sampleID2 . "\t"
			  . $resultsTempArr[2] . "\t"
			  . $resultsTempArr[3] . "\t"
			  . $sampleToGroup{$sampleID1} . "\t"
			  . $sampleToGroup{$sampleID2} . "\t"
			  . $resultsTempArr[4] . "\t"
			  . $resultsTempArr[5] . "\t"
			  . $resultsTempArr[6] . "\n";
		}
	}
}

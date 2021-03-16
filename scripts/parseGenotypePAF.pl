#!/usr/bin/perl

use warnings;
use strict;
use diagnostics;
use IO::File;

#TODO preserve allele order information?

#HERV-C4	6387	0	6387	+	HG01978#2#h2tg000049l:2468298-2626522	158225	52002	58389	6385	6387	0	NM:i:2	ms:i:6367	AS:i:6367	nn:i:0	tp:A:S	cm:i:632	s1:i:6331	de:f:0.0003	rl:i:0	cg:Z:6387M
#alleleID->totalCount
my %alleles;

#sample->alleleIDs->count
my %sampleAlleles;

my $fh      = new IO::File( $ARGV[0], "r" );
my $line    = $fh->getline();
my @tempArr = split( /\t/, $line );

#my @currentSampInfo = split(/#/,$tempArr[5]);
#my $currentID = $tempArr[5];
$alleles{ $tempArr[0] } += 1;
$sampleAlleles{ $tempArr[5] }->{ $tempArr[0] } += 1;

$line = $fh->getline();

while ($line) {
	chomp($line);
	@tempArr = split( /\t/, $line );
	$alleles{ $tempArr[0] } += 1;
	$sampleAlleles{ $tempArr[5] }->{ $tempArr[0] } += 1;
	$line = $fh->getline();
}
$fh->close();

#Print header
print "sample\tcount\ttype\n";
#foreach my $alleleID ( keys %alleles ) {
#	print "\t" . $alleleID;
#}
#print "\n";

#print information
foreach my $sampleID ( keys %sampleAlleles ) {
	foreach my $alleleID ( keys %alleles ) {
		if ( exists( $sampleAlleles{$sampleID}->{$alleleID} ) ) {
			print $sampleID. "\t"
			  . $sampleAlleles{$sampleID}->{$alleleID} ."\t"
			  . $alleleID . "\n";
		}
		else {
			print $sampleID. "\t0\t" . $alleleID . "\n";
		}
	}
}


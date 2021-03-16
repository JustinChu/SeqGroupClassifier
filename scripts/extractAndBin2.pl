#!/usr/bin/perl

use warnings;
use strict;
use diagnostics;
use IO::File;

#CHM13#0#chr6	31786100	31911587	*	0	+
my $fh   = new IO::File( $ARGV[0], "r" );
my $line = $fh->getline();

my %lengths;

while ($line) {
	chomp($line);
	my @tempArr1  = split( /\t/, $line );
	my @tempArr2  = split( /#/, $tempArr1[0] );
	my $filename = $tempArr2[0];
	if ( $tempArr2[1] eq '1' ) {
		$filename .= '.paternal';
	}
	elsif ( $tempArr2[1] eq '2' ) {
		$filename .= '.maternal';
	}
	$filename .= ".fa.gz";
	
	my $region = $tempArr1[0] . ":" . ($tempArr1[1]+1) . "-" . $tempArr1[2];

	#CHM13#0#chr6	31786100	31911587	*	0	+
	#HG00438#2#h2tg000042l	24395718	24488464	*	0	-
	my $cmd =
		'samtools faidx '
	  . $ARGV[1] . '/'. $filename . ' '
	  . $region . ' > '
	  . 'temp.fa';
	#extract region
	system($cmd);
	
	my $length = int(($tempArr1[2] - $tempArr1[1] - 1)/1000);
	
	#bin sequences by length
	$lengths{$length} .= `cat temp.fa`;
	system("rm temp.fa");

	$line = $fh->getline();
}
$fh->close();

foreach my $key (keys %lengths){
	my $fhOut   = new IO::File( $key . "k.fa", "w" );
	$fhOut->write($lengths{$key});
	$fhOut->close();
}

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
	my @tempArr  = split( /#/, $line );
	my $filename = $tempArr[0];
	if ( $tempArr[1] eq '1' ) {
		$filename .= '.paternal';
	}
	elsif ( $tempArr[1] eq '2' ) {
		$filename .= '.maternal';
	}
	$filename .= ".fa.gz";

	#HG00438#2#h2tg000042l:24395719-24488464
	my $cmd =
		'samtools faidx '
	  . $ARGV[1] . '/'. $filename . ' '
	  . $line . ' > '
	  . 'temp.fa';
	#extract region
	system($cmd);
	
	my @tempArr2  = split( /:|-/, $line );
	my $length = int(($tempArr2[2] - $tempArr2[1])/1000);
	
	#bin sequences by length
	$lengths{$length} .= `cat temp.fa`;
	system("rm temp.fa");

	$line = $fh->getline();
}
$fh->close();

foreach my $key (keys %lengths){
	my $fhOut   = new IO::File( $key . "k.fa", "w" );
	$fhOut->write($lengths{$key});
	$fhOut->close();;
}

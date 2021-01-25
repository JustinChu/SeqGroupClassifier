#!/usr/bin/perl

#given reference genomeID and file -> extract from master bed file

use warnings;
use strict;
use diagnostics;
use IO::File;

#CHM13#0#chr6	31786100	31911587	*	0	+
my $fh   = new IO::File( $ARGV[0], "r" );
my $line = $fh->getline();

while ($line) {
	chomp($line);
	my @tempArr = split( /#/, $line );
#	my $tmpFH   = new IO::File( $1 , "w" );
#	$tmpFH->write($)
	my $cmd = 'echo "' . $tempArr[2] . '" >> ' . $tempArr[0] . '.bed';
	system($cmd);
	$line = $fh->getline();
}
$fh->close();
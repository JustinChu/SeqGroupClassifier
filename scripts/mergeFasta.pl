#!/usr/bin/perl

#given reference genomeID and file -> extract from master bed file

use warnings;
use strict;
use diagnostics;
use IO::File;

#CHM13#0#chr6	31786100	31911587	*	0	+
my $fh   = new IO::File( $ARGV[0], "r" );
my $line = $fh->getline();

my %fastaSeq;

while ($line) {
	chomp($line);
	my @tempArr = split( /#/, $line ); 	
	my $newID = $tempArr[0] . "#" . $tempArr[1]. "#";
	
	$line = $fh->getline();
	my $fastaSeq = "";
	
#	if($line =~ /^>/ ){
#		$line = $fh->getline();
#		next;
#	}
	
	while ($line && $line !~ /^>/ ) {
		chomp($line);
		$fastaSeq .= $line;
		$line = $fh->getline();
	}
	if ( exists $fastaSeq{$newID} ) {
		$fastaSeq{$newID} .= "N" . $fastaSeq;
	}
	else {
		$fastaSeq{$newID} = $fastaSeq;
	}
}
$fh->close();

foreach my $id ( keys(%fastaSeq) ) {
	print $id . "\n" . $fastaSeq{$id} . "\n";
}

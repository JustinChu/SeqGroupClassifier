#!/usr/bin/perl

#given reference genomeID and file -> extract from master bed file

use warnings;
use strict;
use diagnostics;
use IO::File;

#chr6	31984664	32023790	>s60685	>s60692	>s60686>s60687>s60688>s60689>s60690>s60691:39126:+:GRCh38#0#chr6:31984661:32023790

#CHM13#0#chr6	31786100	31911587	*	0	+
my $fh   = new IO::File( $ARGV[0], "r" );
my $line = $fh->getline();

while ($line) {
	chomp($line);
	if($line =~ /:([^:]+):([^:]+):(\d+):(\d+)?/){
		my $direction = $1;
		my $chr = $2;
		my $start = $3;
		my $end = $4;
		print $chr . "\t" . $start . "\t" . $end  . "\t*\t\t0\t" . $direction . "\n";
	}
	elsif($line =~ /\.?/){{
		
	}
	else{
		die "failed to parse: $line\n";
	}
	$line = $fh->getline();
}
$fh->close();

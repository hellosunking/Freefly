#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <flash.tab.out> <out.prefix> [motif.kmer=4] [min.qual=30]\n\n";
	exit 2;
}

my $kmer = $ARGV[2] || 4;
if( $kmer <= 0 ) {
	print STDERR "ERROR: kmer MUST be a positive integer!\n";
	exit 1;
}
my $minQ = $ARGV[3] || 30;
$minQ += 33;	## Phred 33

my %size;
my %motif;
my ($all, $long, $valid) = ( 0, 0, 0 );

open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	my @l = split /\t/;	## sid seq1 qual1 [seq2 qual2]
	++ $all;

	if( $#l == 2 ) {	## combined
		++ $size{ length($l[1]) };
	} else {	## uncombined
		++ $long;
	}

	my $pass = 1;
	foreach my $i ( 0..$kmer ) {
		if( ord( substr($l[2], $i, 1) ) < $minQ ) {
			$pass = 0;
			last;
		}
	}
	if( $pass ) {
		my $m = substr( $l[1], 0, $kmer );
		++ $valid;
		++ $motif{$m};
	}
}
close IN;

open SIZE, ">$ARGV[1].size" or die "$!";
foreach my $s ( sort {$a<=>$b} keys %size ) {
	print SIZE join("\t", $s, $size{$s}, $size{$s}/$all*100), "\n";
}
print SIZE join("\t", "#Longer", $long, $long/$all*100), "\n";
close SIZE;

open MOTIF, ">$ARGV[1].motif" or die "$!";
my $entropy = 0;
foreach my $m ( sort keys %motif ) {
	my $freq = $motif{$m}/$valid;
	print MOTIF join("\t", $m, $motif{$m}, $freq*100), "\n";
	if( $freq > 1e-10 ) {
		$entropy += $freq * log($freq);
	}
}
print MOTIF join("\t", "#Invalid", $all-$valid, (1-$valid/$all)*100), "\n";
print MOTIF join("\t", "#Diversity", $valid, -$entropy/log(256)), "\n";
close MOTIF;


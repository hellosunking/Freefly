#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL
#

use strict;
use warnings;

if( $#ARGV < 2 ) {
	print STDERR "\nUsage: $0 <genome.fasta> <in.sam> <output.prefix> [kmer=4] [min.mapQ=30] [autoSome.only=1]\n";
	print STDERR "\nCalculate size and 5'-kmer end motifs.\n\n";
	exit 2;
}

my $kmer = $ARGV[3] || 4;
my $min_mapQ = $ARGV[4] || 30;
my $autosome = (defined $ARGV[5]) ? $ARGV[5] : 1;

my $ACGT_only = 1;	## only consider A/C/G/T, discard all other letters (N,Y,W,etc.)

my $g = load_genome( $ARGV[0] );

if( $ARGV[1] =~ /\.bam$/ ) {
	open IN, "samtools view $ARGV[1] |" or die( "$!" );
} else {
	open IN, "$ARGV[1]" or die( "$!" );
}

my (%motif, %size);
my ($valid_motif, $valid_size) = (0, 0);
while( <IN> ) {
	chomp;
	my @l = split /\t/;	##V06600080699    163     chr10   60818   255     95M     =       60949   225     AACCAACCCAAATGT
	next unless exists $g->{$l[2]} && $l[1] & 2;	## chr exists and the fragment is properly aligned
	next if ($l[4]<$min_mapQ) || ($l[1] & 0xF00);	## low mapQ, 2nd alignment, QC failed, or PCR duplicate
	next if $l[1] & 16;								## reversed strand
	next if $autosome && $l[2]!~/^chr\d+$/;

	## size
	++ $size{ $l[8] };
	++ $valid_size;

	## motif
	my $m = substr( $g->{$l[2]}, $l[3], $kmer );
	if( (!$ACGT_only) || $m=~/^[ACGT]+$/ ) {
		$motif{$m} ++;
		++ $valid_motif;
	}
}
close IN;

open OUT, ">$ARGV[2].size" or die( "$!" );
print OUT "#size\ttCount\tFrequency%\n";
foreach my $i ( sort {$a<=>$b} keys %size ) {
	print OUT join("\t", $i, $size{$i}, $size{$i}/$valid_size*100), "\n";
}
close OUT;

open OUT, ">$ARGV[2].motif" or die( "$!" );
my $entropy = 0;
print OUT "#Motif\tCount\tFrequency%\n";
foreach my $i ( sort keys %motif ) {
	my $freq = $motif{$i}/$valid_motif;
	print OUT join("\t", $i, $motif{$i}, $freq*100), "\n";
	if( $freq > 1e-10 ) {
		$entropy += $freq * log($freq);
	}
}
print OUT join("\t", "#Normalized.entropy",$valid_motif , -$entropy/log(256)), "\n";
close OUT;

sub load_genome {
	my $fasta = shift;

	print STDERR "Loading genome $fasta ...\n";
	my %g;
	my $chr = 'NULL';
	if( $fasta =~ /\.gz$/ ) {
		open IN, "zcat $fasta |" or die("$!");
	} else {
		open IN, "$fasta" or die("$!");
	}
	while( <IN> ) {
		chomp;
		if( /^>(\S+)/ ) {
			$chr = $1;
			$g{$chr} = "X";	## placeholder, as sam format uses 1-base
		} else {
			$g{$chr} .= uc $_;
		}
	}
	close IN;
	return \%g;
}


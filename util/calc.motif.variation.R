#
# Author: Kun Sun @ SZBL
#
# R script for 
#

argv = commandArgs(T);
if( length(argv) != 1 ) {
	print( 'usage: R --slave --args <in.motif> < plot.R' );
	q();
}

motif = read.table( argv[1], row.names=1 );
#Motif	#count	%motif
#AAAA	513466	0.6562
motif = motif[ grepl( "^[ACGT]+$", row.names(motif)), ];

m = motif[,2] / 100;	## extended mode
m.zero = which( m < 1e-16 );
if( length(m.zero) != 0  ) {
	m = m[ - m.zero ];
}
print( paste(argv[1], -sum(m*log2(m))/8, sep=" ") )


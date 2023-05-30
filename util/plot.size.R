#
# Author: ahfyth
#
# R script for 
#

argv = commandArgs(T);
if( length(argv) < 1 ) {
	print( 'usage: R --slave --args <sid> < plot.R' );
	q();
}

sid = argv[1];
pdf( paste0(sid, ".size.pdf") );

b2r_colors = colorRampPalette(c("blue", "red"))(10);
rand_color = b2r_colors[ sample.int(10,1) ];

size = read.table( paste0(sid, ".size") );
plot( size[,3] ~ size[,1], type='l', lwd=2, col=rand_color,
		xlab="Fragment size (bp)", ylab="Frequency (%)", cex.lab=1.2,
		main=paste("Size distribution for ", argv[1], sep=""),
		col.main=rand_color );


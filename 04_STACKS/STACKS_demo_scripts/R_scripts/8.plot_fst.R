#!/usr/bin/env Rscript

sigma='100kb'

x = read.delim('populations.fst_shallow-deep.tsv')

# The large chromosomes.
 gac_groups = "un"
 # c("groupI","groupII","groupIII","groupIV","groupV","groupVI","groupVII","groupVIII","groupIX","groupX","groupXI","groupXII","groupXIII","groupXIV","groupXV","groupXVI","groupXVII","groupXVIII","groupXIX","groupXX","groupXXI")

# Print FST, one page per chromosome.
pdf('fst.pdf')
 	xmax_mb = max(x$BP) %/% 1e6 + 1

	plot(NULL,
		xlim=c(0, xmax_mb*1e6),
		ylim=c(0, max(x$AMOVA.Fst)),
		xlab='Position (Mb)',
		ylab='Fst',
		xaxt='n', # Don't draw the x axis right away.
		main=paste('Fst, deep (SLucia+POS) vs. freshwater (BORD+FAV+PTF+SAV)\n', sep='')#, sigma=', sigma)
		)
	axis(1, at=0:xmax_mb*1e6, labels=0:xmax_mb)

	# Add the SNP-wise FST values.
 	points(x$BP, x$AMOVA.Fst, pch=3, col='grey50', cex=.6)
	
	# Add the smoothed FST.
 	# lines(x$BP, x$Smoothed.AMOVA.Fst, lwd=2, col='blue')
null=dev.off()

#!/usr/bin/env Rscript

library(adegenet)

# For more information, see esp. Jombard's 'An introduction to adegenet 2.0.0',
# Chapter 6 'Multivariate analysis'.

# Load the genotypes.
x = read.genepop('./batch_1.gen')

# Load the popmap and update the genepop object.
popmap = read.table('./popmap.tsv', col.names=c('sample', 'pop'), row.names='sample')
pop(x) = popmap[indNames(x), 'pop']

# Do the PCA.
x.scaled = scaleGen(x, NA.method='mean')
x.pca = dudi.pca(x.scaled, scannf=F, nf=10)

# Plot the PCA.
pdf('pca.pdf')
s.class(x.pca$li, xax=1, yax=2, pop(x), col=rainbow(nPop(x)))
add.scatter.eig(x.pca$eig[1:10], pos="topleft", xax=1, yax=2)
null=dev.off()

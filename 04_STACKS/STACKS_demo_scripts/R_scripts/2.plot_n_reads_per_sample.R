#!/usr/bin/env Rscript

d = read.delim('./n_reads_per_sample.tsv')

pdf('./n_reads_per_sample.pdf', height=12, width=12)
dotchart(sort(setNames(d$n_reads, d$X.sample)), xlim=c(0, max(d$n_reads)))
hist(d$n_reads, nclass=20)

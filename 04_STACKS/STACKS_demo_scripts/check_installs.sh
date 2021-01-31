#!/bin/bash

echo '# Perl:'
command -v perl || perl

echo '# R:'
command -v R || R

echo '# Samtools:'
command -v samtools || samtools

echo '# BWA:'
command -v bwa || bwa

echo '# Stacks:'
command -v cstacks || cstacks

echo '# ADEgenet (trying to load the package):'
R --slave -e "library(adegenet);"

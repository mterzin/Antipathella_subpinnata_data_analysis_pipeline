#!/bin/bash
top=$(readlink -f $(dirname $0)/..)

# STEP 12: Download the genome
# ==========
echo "Downloading the genome..."
cd $top/genome
wget "ftp://ftp.ensembl.org/pub/release-87/fasta/gasterosteus_aculeatus/dna/Gasterosteus_aculeatus.BROADS1.dna.toplevel.fa.gz" || exit

# STEP 13: Create the BWA database
# ==========
echo "Creating BWA database..."
genome_fa=Gasterosteus_aculeatus.BROADS1.dna.toplevel.fa.gz
mkdir bwa
bwa index -p bwa/gac $genome_fa &> bwa/bwa_index.oe

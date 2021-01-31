#!/bin/bash
top=$(readlink -f $(dirname $0)/..)
echo "$top"

	######################################
	# Set your filtering parameters here #
	######################################
	
	# Minimum stack depth/coverage
	M=3
	# -p : minimum number of populations a locus must be present in to process a locus
	p=1
	# -R : minimum percentage of individuals in a metapopulation required to process a locus for that population.
	R=0.8
	# -r : minimum percentage of individuals in a population required to process a locus for that population.
	# --max_obs_het : maximum observed heterozygosity required to process a nucleotide site at a locus.
	# --min-maf : minimum minor allele frequency required to process a nucleotide site at a locus (0 < min_maf < 0.5)
	maf=0.1
	# --min-mac : minimum minor allele count required to process a nucleotide site at a locus (integer).
	mac=2

#############################################################################################
#  I don't think I can plot Fst with this script, because I don't have a reference genome   #
# I will erase steps 23 and 24 because of that, but I can still see it from the demo script #
#############################################################################################

# ==========

# STEP 22: Create a new "deep vs shallow" popmap.
cd $top/info
sed -r 's/\t(POS|SLucia)/\tdeep/; s/\t(BORD|FAV|PTF|SAV)/\tshallow/;' popmap.tsv \
	> popmap.deep_shallow.tsv
	# I only kept this because the bash command is useful, I don't think I will be using this file

# PCA
# ==========

cd $top/stacks.denovo/stacks.M$M.Mediterranean_POSADA_concatenated

# STEP 25: Filter the genotypes and export them in GENEPOP format (with populations).
echo "Running populations with --genepop..."
popmap=../../info/popmap_no_samples_with_low_reads_POSADA_concatenated.tsv
mkdir pops.p$p-R$R-maf$maf-mac$mac-POSADA_concatenated
out_dir=pops.p$p-R$R-maf$maf-mac$mac-POSADA_concatenated
log_file=pops.p$p-R$R-maf$maf-mac$mac-POSADA_concatenated.oe
export PATH=/usr/arb/bin:/home/usr7eco1/Programs/SHRiMP/bin:/home/usr7eco1/Programs/bin:/home/usr7eco1/Paolo/Bacteria_exp_field_copy/Analyses/ncbi-blast-2.9.0+/bin:/home/usr7eco1/Programs/Stacks/bin:/home/usr7eco1/perl5/bin:/usr/arb/bin:/home/usr7eco1/Programs/SHRiMP/bin:/home/usr7eco1/Programs/bin:/home/usr7eco1/Paolo/Bacteria_exp_field_copy/Analyses/ncbi-blast-2.9.0+/bin:/home/usr7eco1/Programs/Stacks/bin:/home/usr7eco1/perl5/bin:/usr/lib64/qt-3.3/bin:/usr/lib64/ccache:/usr/local/bin:/usr/local/sbin:/usr/bin:/usr/sbin:/bin:/sbin:/home/usr7eco1/.local/bin:/home/usr7eco1/bin${PATH} && /home/usr7eco1/Programs/stacks-2.3e/populations -P ./ -M $popmap -O $out_dir \
	-p $p -R $R --min_mac $mac --min_maf $maf --write_single_snp --hwe --fstats \
	--genepop --vcf --structure \
	&> $log_file
	
# STEP 26: Copy, move, or "symlink" the GENEPOP file (ADEgenet expects a `.gen` suffix).
# We also symlink the population map so that we can find it easily in the R script.
cd $top/stacks.denovo/stacks.M$M.Mediterranean_POSADA_concatenated/$out_dir
mkdir adegenet
cd adegenet
ln -s ../populations.snps.genepop batch_1.gen
popmap=../../$popmap
ln -s $popmap .

##################################################################################################################
# IMPORTANT
# Since there is no adegenet package on the cluster I stop here and run the analysis in R on my laptop
##################################################################################################################

# STEP 27, STEP28: Perform the PCA (with R).
# echo "Performing and plotting the PCA using ADEgenet..."
# $top/demo_scripts/R_scripts/8.plot_pca.R

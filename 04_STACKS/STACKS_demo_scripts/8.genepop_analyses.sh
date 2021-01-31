#!/bin/bash
top=$(readlink -f $(dirname $0)/..)

# Plotting FST along one chromosome.
# ==========

# STEP 22: Create a new "oceanic vs freshwater" popmap.
cd $top/info
sed -r 's/\t(cs|sj)/\toceanic/; s/\t(pcr|wc)/\tfreshwater/;' popmap.tsv \
	> popmap.oceanic_freshwater.tsv

# STEP 23: Run populations.
cd $top/stacks.ref/rxstacks
echo "Running populations with --fstats -k..."
popmap=../../info/popmap.oceanic_freshwater.tsv
out_dir=pops.oceanic_freshwater
log_file=pops.oceanic_freshwater.oe
populations -P ./ -M $popmap -O $out_dir \
	-p 2 -r 0.80 --min_maf 0.05 --max_obs_het 0.70 \
    --fstats -k --sigma 100000 \
    &> $log_file

# STEP 24: Plot the results with R.
echo "Plotting the FST scan..."
cd $out_dir
$top/demo_scripts/R_scripts/8.plot_fst.R

# PCA
# ==========

cd $top/stacks.denovo/rxstacks

# STEP 25: Filter the genotypes and export them in GENEPOP format (with populations).
echo "Running populations with --genepop..."
popmap=../../info/popmap.tsv
out_dir=pops.r80-maf05-het70
log_file=pops.r80-maf05-het70.oe
populations -P ./ -M $popmap -O $out_dir \
	-p 4 -r 0.80 --min_maf 0.05 --max_obs_het 0.70 \
	--genepop \
	&> $log_file

# STEP 26: Copy, move, or "symlink" the GENEPOP file (ADEgenet expects a `.gen` suffix).
# We also symlink the population map so that we can find it easily in the R script.
cd $out_dir
mkdir adegenet
cd adegenet
ln -s ../batch_1.genepop batch_1.gen
popmap=../../$popmap
ln -s $popmap .

# STEP 27, STEP28: Perform the PCA (with R).
echo "Performing and plotting the PCA using ADEgenet..."
$top/demo_scripts/R_scripts/8.plot_pca.R

#!/bin/bash
top=$(readlink -f $(dirname $0)/..)

# STEP 16, STEP 17-A-i: Change directory.
cd $top/stacks.denovo

# STEP 17-A-ii: Run ustacks on every sample.
M=4
n=4
index=1
for sample in $(cut -f1 ../info/popmap.tsv) ;do
	fq_file=../cleaned/$sample.fq.gz
	logfile=$sample.ustacks.oe
	ustacks -f $fq_file -i $index -o ./ -M $M -m 3 &> $logfile
	index=$(( $index + 1 ))
done

# STEP 17-A-iii: Check that all ustacks runs have completed.
echo "Checking that all ustacks runs have completed..."
ls *.ustacks.oe | wc
wc -l *.ustacks.oe
grep -iE '\b(err|e:|warn|w:|fail|abort)' *.ustacks.oe
grep -L 'ustacks is done' *.ustacks.oe

# STEP 17-A-iv: Extract the sample coverages from the `*.ustacks.oe` files
# and check the values are consistent with previously obtained coverages.

# STEP 17-A-v: Pick the samples to include in the catalog.
# For each population, we pick the 10 samples with the highest coverage (40
# samples total), excluding the high-coverage outlier sj_1484.07.
n_per_pop=10
excluded_samples=sj_1484.07
for pop in $(cut -f2 ../info/popmap.tsv | sort -u) ;do
	sort -k2,2nr ../info/n_reads_per_sample.tsv \
		| grep -v $excluded_samples \
		| grep "^$pop" \
		| head -n $n_per_pop \
		| cut -f1 \
		| sed -r 's/([A-Za-z]+)(.+)/\1\2\t\1/'
done > ../info/popmap.catalog.tsv

# STEP 17-A-vi: Run cstacks
cstacks -P ./ -M ../info/popmap.catalog.tsv -n $n &> cstacks.oe

# STEP 17-A-vii: Run sstacks on every sample
for sample in $(cut -f1 ../info/popmap.tsv) ;do
	logfile=$sample.sstacks.oe
	sstacks -c ./batch_1 -s ./$sample -o ./ &> $logfile
done

# (Check that all sstacks runs have completed.)

# STEP 18: We use rxstacks to improve calls. We create a new subdirectory.
mkdir rxstacks
cd rxstacks

# STEP 19: Run rxstacks.
rxstacks -P ../ -o ./ --prune_haplo --conf_lim=0.10 &> rxstacks.oe

# STEP 20: Re-run cstacks and sstacks
# (Same commands as above, with updated paths.)
cstacks -P ./ -M ../../info/popmap.catalog.tsv -n $n &> cstacks.oe
for sample in $(cut -f1 ../../info/popmap.tsv) ;do
	logfile=$sample.sstacks.oe
	sstacks -c ./batch_1 -s ./$sample -o ./ &> $logfile
done

# STEP 21: Filter genotypes with populations.
min_samples=0.80
min_maf=0.05
max_obs_het=0.70
populations -P ./ -r $min_samples --min_maf=$min_maf --max_obs_het=$max_obs_het &> populations.oe

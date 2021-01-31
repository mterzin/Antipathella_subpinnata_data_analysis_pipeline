#!/bin/bash
top=$(readlink -f $(dirname $0)/..)

# STEP 14: Pick twelve representative samples.
echo -n "\
cs_1335.01	cs
cs_1335.13	cs
cs_1335.15	cs
pcr_1193.10	pcr
pcr_1193.11	pcr
pcr_1211.04	pcr
sj_1483.06	sj
sj_1484.07	sj
sj_1819.36	sj
wc_1218.04	wc
wc_1221.01	wc
wc_1222.02	wc
" > $top/info/popmap.test_samples.tsv

# STEP 15-A-i: Chose parameter combinations to survey.
M_values="1 2 3 4 5 6 7 8 9"

# STEP 15-A-ii: Change directory.
cd $top/tests.denovo

# STEP 15-A-iii: Create subdirectories.
for M in $M_values ;do
	mkdir stacks.M$M
done

# STEP 15-A-iv: Run denovo_map on the subset of samples.
popmap=../info/popmap.test_samples.tsv
for M in $M_values ;do
	n=$M
	m=3
	echo "Running Stacks for M=$M, n=$n..."
	reads_dir=../cleaned
	out_dir=stacks.M$M
	log_file=$out_dir/denovo_map.oe
	denovo_map.pl --samples $reads_dir -O $popmap -o $out_dir -b 1 -S -M $M -n $n -m $m &> $log_file
done

# STEP 15-A-v: Check that all runs have completed.
echo "Checking that all denovo_map runs have completed..."
ls stacks.M*/denovo_map.oe | wc
wc -l stacks.M*/denovo_map.oe
grep -iE '\b(err|e:|warn|w:|fail|abort)' stacks.M*/denovo_map.oe stacks.M*/denovo_map.log
grep -L 'denovo_map\.pl is done' stacks.M*/denovo_map.log

# STEP 15-A-vi: Check coverages (in file `stacks.M1/denovo_map.log`).

# STEP 15-A-vii: Run populations with '-r 0.80' (loci present in 80% of samples)
for M in $M_values ;do
	stacks_dir=stacks.M$M
	out_dir=$stacks_dir/populations.r80
	mkdir $out_dir
	log_file=$out_dir/populations.oe
	populations -P $stacks_dir -O $out_dir -r 0.80 &> $log_file
done

# STEP 15-A-viii: Compare the results obtained with different parameters
mkdir results
cd results

# Extract the SNPs-per-locus distributions. These distributions
# are reported in the log of populations.
echo "Tallying the numbers..."
echo -e '#par_set\tM\tn\tm\tn_snps\tn_loci' > n_snps_per_locus.tsv
for M in $M_values ;do
	n=$M
	m=3
	log_file=../stacks.M$M/populations.r80/batch_1.populations.log

	# Extract the numbers for this parameter combination.
	sed -n '/^#n_snps\tn_loci/,/^[^0-9]/ p' $log_file | grep -E '^[0-9]' > $log_file.snps_per_loc

	# Cat the content of this file, prefixing each line with information on this
	# parameter combination.
	line_prefix="M$M-n$n-m$m\t$M\t$n\t$m\t"
	sed -r "s/^/$line_prefix/" $log_file.snps_per_loc >> n_snps_per_locus.tsv
done

# Plot the results with R.
echo "Plotting the number of loci..."
$top/demo_scripts/R_scripts/4.plot_n_loci.R
echo "Plotting the distribution of the number of SNPs..."
$top/demo_scripts/R_scripts/4.plot_n_snps_per_locus.R

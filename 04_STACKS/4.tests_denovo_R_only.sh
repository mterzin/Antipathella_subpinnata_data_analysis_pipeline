#!/bin/bash
top=$(readlink -f $(dirname $0)/..)
echo "$top"

# STEP 14: Pick 15 representative samples.
echo -n "\
T_AsBORD_1_S12	BORD
T_AsBORD_2_S14	BORD
T_AsBORD_3_S17	BORD
T_AsBORD_4_S20	BORD
T_AsBORD_5_S21	BORD
T_AsBORD_6_S25	BORD
T_AsBORD_7_S26	BORD
T_AsBORD_8_S27	BORD
T_AsBORD_9_S2	BORD
T_AsBORD_10_S28	BORD
T_AsFAV_1_S5	FAV
T_AsFAV_2_S29	FAV
T_AsFAV_3_S6	FAV
T_AsFAV_4_S1	FAV
T_AsFAV_5_S7	FAV
T_AsFAV_6_S30	FAV
T_AsFav_7_S32	FAV
T_AsFAV_8_S33	FAV
T_AsFAV_9_S34	FAV
T_AsFAV_10_S35	FAV
T_AsFAV_11_S36	FAV
T_AsPOS_1_S73	POS
T_AsPOS_2r_S79	POS
T_AsPOS_3r_S80	POS
T_AsPOS_4r_S85	POS
T_AsPOS_5_S77	POS
T_AsPTF_1_S44	PTF
T_AsPTF_2_S45	PTF
T_AsPTF_3_S60	PTF
T_AsPTF_4_S8	PTF
T_AsPTF_5_S52	PTF
T_AsPTF_6_S71	PTF
T_AsPTF_7_S61	PTF
T_AsPTF_8_S53	PTF
T_AsPTF_9_S66	PTF
T_AsPTF_10_S92	PTF
T_AsSAV_1_S82	SAV
T_AsSAV_2_S69	SAV
T_AsSAV_3_S84	SAV
T_AsSAV_4_S93	SAV
T_AsSAV_6_S87	SAV
T_AsSAV_7_S94	SAV
T_AsSAV_8r_S91	SAV
T_AsSAV_9_S68	SAV
T_AsSAV_10_S64	SAV
T_AsSLucia_1_S63	SLucia
T_AsSLucia_2r_S65	SLucia
T_AsSLucia_3A_S62	SLucia
T_AsSLucia_5_S67	SLucia
" > $top/info/popmap.tsv

# STEP 15-A-i: Chose parameter combinations to survey.
#M_values="1 2 3 4 5 6 7 8 9"
 
# STEP 15-A-ii: Change directory.
cd $top/tests.denovo

# STEP 15-A-iii: Create subdirectories.
#for M in $M_values ;do
#	mkdir stacks.M$M
#done

# STEP 15-A-iv: Run denovo_map on the subset of samples.
popmap=$top/info/popmap.tsv
echo "File containing all the samples to do a test run is: $popmap"

#for M in $M_values ;do
#	n=$M
#	m=3
#	echo "Running Stacks for M=$M, n=$n..."
#	reads_dir=$top/cleaned
#	echo "Input files (quality filtered) are in $reads_dir"
#	out_dir=$top/tests.denovo/stacks.M$M
#	echo "This is the output directory: $out_dir"
#	log_file=$out_dir/denovo_map.oe
#	echo "Log file is $log_file"
#	export PATH=/usr/arb/bin:/home/usr7eco1/Programs/SHRiMP/bin:/home/usr7eco1/Programs/bin:/home/usr7eco1/Paolo/Bacteria_exp_field_copy/Analyses/ncbi-blast-2.9.0+/bin:/home/usr7eco1/Programs/Stacks/bin:/home/usr7eco1/perl5/bin:/usr/arb/bin:/home/usr7eco1/Programs/SHRiMP/bin:/home/usr7eco1/Programs/bin:/home/usr7eco1/Paolo/Bacteria_exp_field_copy/Analyses/ncbi-blast-2.9.0+/bin:/home/usr7eco1/Programs/Stacks/bin:/home/usr7eco1/perl5/bin:/usr/lib64/qt-3.3/bin:/usr/lib64/ccache:/usr/local/bin:/usr/local/sbin:/usr/bin:/usr/sbin:/bin:/sbin:/home/usr7eco1/.local/bin:/home/usr7eco1/bin${PATH} && /home/usr7eco1/Programs/stacks-2.3e/scripts/denovo_map.pl --samples $reads_dir --popmap $popmap -o $out_dir -M $M -n $n -m $m &> $log_file
#done

# STEP 15-A-v: Check that all runs have completed.
#echo "Checking that all denovo_map runs have completed..."
#ls stacks.M*/denovo_map.oe | wc
#wc -l stacks.M*/denovo_map.oe
#grep -iE '\b(err|e:|warn|w:|fail|abort)' stacks.M*/denovo_map.oe stacks.M*/denovo_map.log
#grep -L 'denovo_map\.pl is done' stacks.M*/denovo_map.log

# STEP 15-A-vi: Check coverages (in file `stacks.M1/denovo_map.log`).

# STEP 15-A-vii: Run populations with '-r 0.80' (loci present in 80% of samples)
#for M in $M_values ;do
#	stacks_dir=$top/tests.denovo/stacks.M$M
#	echo "This directory contains STACKS files: $stacks_dir"
#	out_dir=$stacks_dir/populations.r80
#	echo "This is the output directory: $out_dir"
#	mkdir $out_dir
#	log_file=$out_dir/populations.oe
#	echo "Log file is $log_file"
#	export PATH=/usr/arb/bin:/home/usr7eco1/Programs/SHRiMP/bin:/home/usr7eco1/Programs/bin:/home/usr7eco1/Paolo/Bacteria_exp_field_copy/Analyses/ncbi-blast-2.9.0+/bin:/home/usr7eco1/Programs/Stacks/bin:/home/usr7eco1/perl5/bin:/usr/arb/bin:/home/usr7eco1/Programs/SHRiMP/bin:/home/usr7eco1/Programs/bin:/home/usr7eco1/Paolo/Bacteria_exp_field_copy/Analyses/ncbi-blast-2.9.0+/bin:/home/usr7eco1/Programs/Stacks/bin:/home/usr7eco1/perl5/bin:/usr/lib64/qt-3.3/bin:/usr/lib64/ccache:/usr/local/bin:/usr/local/sbin:/usr/bin:/usr/sbin:/bin:/sbin:/home/usr7eco1/.local/bin:/home/usr7eco1/bin${PATH} && /home/usr7eco1/Programs/stacks-2.3e/populations -P $stacks_dir -O $out_dir -r 0.80 -M $popmap &> $log_file
#done

# STEP 15-A-viii: Compare the results obtained with different parameters
#mkdir $top/tests.denovo/results
cd $top/tests.denovo/results

# Extract the SNPs-per-locus distributions. These distributions are reported in the log of populations.
#echo "Tallying the numbers..."
#echo -e '#par_set\tM\tn\tm\tn_snps\tn_loci' > n_snps_per_locus.tsv

#for M in $M_values ;do
#	n=$M
#	m=3
#	log_file=$top/tests.denovo/stacks.M$M/populations.r80/populations.log.distribs
#	echo "Log file is $log_file"

	# Extract the numbers for this parameter combination. We are interested in the number of snps and loci after filtering (snps_per_loc_postfilters)
#	export PATH=/usr/arb/bin:/home/usr7eco1/Programs/SHRiMP/bin:/home/usr7eco1/Programs/bin:/home/usr7eco1/Paolo/Bacteria_exp_field_copy/Analyses/ncbi-blast-2.9.0+/bin:/home/usr7eco1/Programs/Stacks/bin:/home/usr7eco1/perl5/bin:/usr/arb/bin:/home/usr7eco1/Programs/SHRiMP/bin:/home/usr7eco1/Programs/bin:/home/usr7eco1/Paolo/Bacteria_exp_field_copy/Analyses/ncbi-blast-2.9.0+/bin:/home/usr7eco1/Programs/Stacks/bin:/home/usr7eco1/perl5/bin:/usr/lib64/qt-3.3/bin:/usr/lib64/ccache:/usr/local/bin:/usr/local/sbin:/usr/bin:/usr/sbin:/bin:/sbin:/home/usr7eco1/.local/bin:/home/usr7eco1/bin${PATH} && /home/usr7eco1/Programs/stacks-2.3e/scripts/stacks-dist-extract $log_file snps_per_loc_postfilters > $log_file.snps_per_loc
	
	############################################################################
	# This was in the original demo script although I didn't use it in the end #
	#################################################################################################
	# sed -n '/^#n_snps\tn_loci/,/^[^0-9]/ p' $log_file | grep -E '^[0-9]' > $log_file.snps_per_loc #
	#################################################################################################

	# Cat the content of this file, prefixing each line with information on this
	# parameter combination.
#	line_prefix="M$M-n$n-m$m\t$M\t$n\t$m\t"
#	sed -r "s/^/$line_prefix/" $log_file.snps_per_loc >> n_snps_per_locus.tsv
	# This file needs to be modified after by hand. You want to remove the hashtags with the comments
#done

# Plot the results with R.
echo "Plotting the number of loci..."
$top/demo_scripts/R_scripts/4.plot_n_loci.R
echo "Plotting the distribution of the number of SNPs..."
$top/demo_scripts/R_scripts/4.plot_n_snps_per_locus.R

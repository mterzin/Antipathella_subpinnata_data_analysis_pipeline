#!/bin/bash
top=$(readlink -f $(dirname $0)/..)

# STEP 7: Change directory.
cd $top/cleaned

# STEP 8: Clean & demultiplex the data.
echo "Cleaning and demultiplexing the reads..."
for lane in lane1 lane2 lane3 ;do
	barcodes=../info/barcodes.$lane.tsv
	process_radtags -p ../raw/$lane/ -o ./ -b $barcodes -e sbfI --inline_null -cqr &> process_radtags.$lane.oe
done

# STEP 9: Check the per-sample coverages.

# Extract the number of reads from the log of process_radtags.
cd $top/info

echo -e '#sample\tn_reads' > n_reads_per_sample.tsv
for lane in lane1 lane2 lane3 ;do
	# Retrieve the part of the log between 'Barcode...' and the next empty line,
	# then discard the first and last lines and keep the 2nd and 6th columns.
	sed -n '/^Barcode\tFilename\t/,/^$/ p' ../cleaned/process_radtags.$lane.log \
		| sed '1 d; $ d' \
		| cut -f2,6 \
		>> n_reads_per_sample.tsv
done

# Plot these numbers.
echo "Plotting per-sample coverage..."
../demo_scripts/R_scripts/2.plot_n_reads_per_sample.R

# Remove bad samples (if any) from the population map.
echo -n "\
sj_1483.05
sj_1819.31
" > discarded_samples

discarded_regex=$(paste -s -d'|' discarded_samples)
mv popmap.tsv complete_popmap.tsv
grep -vE "^($discarded_regex)\>" complete_popmap.tsv > popmap.tsv

# STEP 10: PCR duplicates (not applicable to the demonstration dataset).
# STEP 11: Read pairs (not applicable to the demonstration dataset).

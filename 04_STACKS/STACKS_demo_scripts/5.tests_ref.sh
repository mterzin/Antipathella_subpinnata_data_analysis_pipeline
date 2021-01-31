#!/bin/bash
top=$(readlink -f $(dirname $0)/..)

# STEP 15-B-i: Choose a temptative alignment method.
# For this example we use BWA with default parameters.

# STEP 15-B-ii: Change directory.
cd $top/tests.ref

# STEP 15-B-iii: Create subdirectories.
mkdir alignments.bwa
mkdir stacks.bwa

# STEP 15-B-iv: Align every test sample.
echo "Running BWA..."
for sample in $(cut -f1 ../info/popmap.test_samples.tsv) ;do
	fq_file=../cleaned/$sample.fq.gz
	bam_file=alignments.bwa/$sample.bam
	log_file=alignments.bwa/$sample.oe
	db=../genome/bwa/gac
	{
		bwa mem -M $db $fq_file | samtools view -b > $bam_file
	} &> $log_file
done

# STEP 15-B-v: Check that alignments have completed.
echo "Checking that all BWA runs have completed..."
ls alignments.bwa/*.oe | wc
wc -l alignments.bwa/*.oe
grep -iE '\b(err|e:|warn|w:|fail|abort)' alignments.bwa/*.oe
tail -n5 alignments.bwa/*.oe | head -n100

# STEP 15-B-vi: Change directory.
cd stacks.bwa

# STEP 15-B-vii: Run pstacks on every test sample.
echo "Running pstacks on default BWA alignments..."
popmap=../../info/popmap.test_samples.tsv
index=1
for sample in $(cut -f1 $popmap) ;do
	bam_file=../alignments.bwa/$sample.bam
	log_file=$sample.pstacks.oe
	pstacks -f $bam_file -i $index -o ./ &> $log_file
	index=$(( $index + 1 ))
done

# STEP 15-B-viii: Check that all pstacks runs have completed.
echo "Checking that all pstacks runs have completed..."
ls *.pstacks.oe | wc
wc -l *.pstacks.oe
grep -iE '\b(err|e:|warn|w:|fail|abort)' *.pstacks.oe
grep -L 'pstacks is done' *.pstacks.oe

# STEP 15-B-ix, STEP 15-B-x, STEP 15-B-xi: Assess the results
mkdir ../results
cd ../results
n_reads_per_sample=../../info/n_reads_per_sample.tsv
echo -e '#sample\tn_primary_alns\tused_reads\tclip_discards\tcoverage' > pstacks_results.tsv
for sample in $(cut -f1 $popmap) ;do
	pstacks_log=../stacks.bwa/$sample.pstacks.oe

	# Initial number of reads
	n_reads=$(grep "^$sample" $n_reads_per_sample | cut -f2)

	# Number of primary alignments that were retained.
	regex='^ +Kept [0-9]+ primary alignments$'
	n_primary_alns=$(grep -E "$regex" $pstacks_log | grep -oE '[0-9]+')
	pct_kept=$(echo "scale=1; $n_primary_alns * 100 / $n_reads" | bc)

	# Percentage of alignments discarded because of clipping
	regex='^ +Skipped [0-9]+ \([0-9.]+%\) excessively soft-clipped primary alignments$'
	pct_clipped=$(grep -E "$regex" $pstacks_log | awk '{print $3}' | tr -d '%()')

	# Mean coverage
	regex='^Created [0-9]+ loci; mean coverage is [0-9.]+ \(stdev: [0-9.]+, max: [0-9]+\).$'
	mean_cov=$(grep -E "$regex" $pstacks_log | awk '{print $7}')

	echo "$sample\t$n_primary_alns\t$pct_kept%%\t$pct_clipped%%\t$mean_cov" >> pstacks_results.tsv
done
